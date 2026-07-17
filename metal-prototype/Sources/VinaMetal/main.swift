import Foundation
import Metal

struct GridParams {
    var dimX: UInt32
    var dimY: UInt32
    var dimZ: UInt32
    var mapCount: UInt32
    var atomCount: UInt32
    var poseCount: UInt32
    var spacing: Float
    var slope: Float
    var originX: Float
    var originY: Float
    var originZ: Float
    var curlV: Float
}

struct PackedFloat4 {
    var x: Float
    var y: Float
    var z: Float
    var w: Float
}

struct SearchParams {
    var lanes: UInt32
    var steps: UInt32
    var seed: UInt32
    var localSteps: UInt32
    var torsionCount: UInt32
    var pairCount: UInt32
    var centerX: Float
    var centerY: Float
    var centerZ: Float
    var padding0: Float = 0
    var spanX: Float
    var spanY: Float
    var spanZ: Float
    var translationMutation: Float
    var rotationMutation: Float
    var temperature: Float
    var gradientStep: Float
    var padding1: Float = 0
}

struct ParsedMap {
    let dims: (Int, Int, Int)
    let spacing: Float
    let center: (Float, Float, Float)
    let values: [Float]
}

struct LigandAtom {
    let lineIndex: Int
    let serial: Int
    let position: PackedFloat4
    let adType: String
}

struct TypedLigand {
    let atoms: [PackedFloat4]
    let xsTypes: [UInt32]
    let originalIndices: [Int]
    let bonds: [[Int]]
}

struct TorsionData {
    var parent: UInt32
    var child: UInt32
    var maskLow: UInt32
    var maskHigh: UInt32
}

struct PairData {
    var a: UInt32
    var b: UInt32
    var typeA: UInt32
    var typeB: UInt32
}

func parseMap(_ url: URL) throws -> ParsedMap {
    let lines = try String(contentsOf: url, encoding: .utf8).split(whereSeparator: \.isNewline).map(String.init)
    guard lines.count > 6 else { throw NSError(domain: "VinaMetal", code: 1, userInfo: [NSLocalizedDescriptionKey: "Invalid map: \(url.path)"]) }
    let spacing = Float(lines[3].split(separator: " ")[1])!
    let elements = lines[4].split(separator: " ").dropFirst().map { Int($0)! + 1 }
    let center = lines[5].split(separator: " ").dropFirst().map { Float($0)! }
    let values = lines.dropFirst(6).map { Float($0)! }
    guard values.count == elements[0] * elements[1] * elements[2] else {
        throw NSError(domain: "VinaMetal", code: 2, userInfo: [NSLocalizedDescriptionKey: "Map value count mismatch"])
    }
    return ParsedMap(dims: (elements[0], elements[1], elements[2]), spacing: spacing, center: (center[0], center[1], center[2]), values: values)
}

func parseLigand(_ url: URL) throws -> (lines: [String], atoms: [LigandAtom]) {
    let lines = try String(contentsOf: url, encoding: .utf8).split(whereSeparator: \.isNewline).map(String.init)
    let atoms = lines.enumerated().compactMap { lineIndex, line -> LigandAtom? in
        guard line.hasPrefix("ATOM") || line.hasPrefix("HETATM") else { return nil }
        let fields = line.split(whereSeparator: \.isWhitespace)
        guard fields.count >= 11, line.utf8.count >= 54, let serial = Int(fields[1]) else { return nil }
        let bytes = Array(line.utf8)
        guard let x = Float(String(decoding: bytes[30..<38], as: UTF8.self).trimmingCharacters(in: .whitespaces)),
              let y = Float(String(decoding: bytes[38..<46], as: UTF8.self).trimmingCharacters(in: .whitespaces)),
              let z = Float(String(decoding: bytes[46..<54], as: UTF8.self).trimmingCharacters(in: .whitespaces)) else { return nil }
        return LigandAtom(lineIndex: lineIndex, serial: serial, position: PackedFloat4(x: x, y: y, z: z, w: 0), adType: String(fields.last!))
    }
    return (lines, atoms)
}

func typedHeavyAtoms(_ ligand: [LigandAtom], mapIndices: [String: Int]) throws -> TypedLigand {
    let radii: [String: Float] = [
        "C": 0.77, "A": 0.77, "N": 0.75, "NA": 0.75, "O": 0.73, "OA": 0.73,
        "P": 1.06, "S": 1.02, "SA": 1.02, "H": 0.37, "HD": 0.37,
        "F": 0.71, "I": 1.33, "Cl": 0.99, "Br": 1.14,
    ]
    func isHydrogen(_ type: String) -> Bool { type == "H" || type == "HD" }
    func isHetero(_ type: String) -> Bool { !["A", "C", "H", "HD"].contains(type) }
    var bonds = Array(repeating: [Int](), count: ligand.count)
    for i in ligand.indices {
        for j in ligand.indices where j > i {
            guard let ri = radii[ligand[i].adType], let rj = radii[ligand[j].adType] else { continue }
            let dx = ligand[i].position.x - ligand[j].position.x
            let dy = ligand[i].position.y - ligand[j].position.y
            let dz = ligand[i].position.z - ligand[j].position.z
            if dx * dx + dy * dy + dz * dz < pow(1.1 * (ri + rj), 2) {
                bonds[i].append(j)
                bonds[j].append(i)
            }
        }
    }
    func mapName(for index: Int) -> String? {
        let type = ligand[index].adType
        if isHydrogen(type) { return nil }
        switch type {
        case "C", "A": return bonds[index].contains { isHetero(ligand[$0].adType) } ? "C_P" : "C_H"
        case "N", "NA":
            let donor = bonds[index].contains { ligand[$0].adType == "HD" }
            let acceptor = type == "NA"
            return acceptor ? (donor ? "N_DA" : "N_A") : (donor ? "N_D" : "N_P")
        case "O", "OA":
            let donor = bonds[index].contains { ligand[$0].adType == "HD" }
            let acceptor = type == "OA"
            return acceptor ? (donor ? "O_DA" : "O_A") : (donor ? "O_D" : "O_P")
        case "S", "SA": return "S_P"
        case "P": return "P_P"
        case "F": return "F_H"
        case "Cl": return "Cl_H"
        case "Br": return "Br_H"
        case "I": return "I_H"
        default: return nil
        }
    }
    let xsIndices: [String: UInt32] = [
        "C_H": 0, "C_P": 1, "N_P": 2, "N_D": 3, "N_A": 4, "N_DA": 5,
        "O_P": 6, "O_D": 7, "O_A": 8, "O_DA": 9, "S_P": 10, "P_P": 11,
        "F_H": 12, "Cl_H": 13, "Br_H": 14, "I_H": 15,
    ]
    var atoms: [PackedFloat4] = []
    var types: [UInt32] = []
    var originalIndices: [Int] = []
    for index in ligand.indices {
        guard let name = mapName(for: index) else { continue }
        guard let mapIndex = mapIndices[name] else {
            throw NSError(domain: "VinaMetal", code: 3, userInfo: [NSLocalizedDescriptionKey: "Missing affinity map \(name) for atom \(ligand[index].serial)"])
        }
        let p = ligand[index].position
        atoms.append(PackedFloat4(x: p.x, y: p.y, z: p.z, w: Float(mapIndex)))
        types.append(xsIndices[name]!)
        originalIndices.append(index)
    }
    return TypedLigand(atoms: atoms, xsTypes: types, originalIndices: originalIndices, bonds: bonds)
}

func parseTorsions(lines: [String], ligand: [LigandAtom], heavyOriginalIndices: [Int]) -> [TorsionData] {
    struct Branch { var parent: Int; var child: Int; var descendants: Set<Int> }
    var branches: [Branch] = []
    var stack: [Int] = []
    for line in lines {
        let fields = line.split(whereSeparator: \.isWhitespace)
        if fields.first == "BRANCH", fields.count >= 3, let parent = Int(fields[1]), let child = Int(fields[2]) {
            branches.append(Branch(parent: parent, child: child, descendants: []))
            stack.append(branches.count - 1)
        } else if fields.first == "ENDBRANCH" {
            _ = stack.popLast()
        } else if (fields.first == "ATOM" || fields.first == "HETATM"), fields.count > 1, let serial = Int(fields[1]) {
            for branchIndex in stack { branches[branchIndex].descendants.insert(serial) }
        }
    }
    let serialToOriginal = Dictionary(uniqueKeysWithValues: ligand.indices.map { (ligand[$0].serial, $0) })
    let originalToHeavy = Dictionary(uniqueKeysWithValues: heavyOriginalIndices.enumerated().map { ($0.element, $0.offset) })
    return branches.compactMap { branch in
        guard let parentOriginal = serialToOriginal[branch.parent], let childOriginal = serialToOriginal[branch.child],
              let parent = originalToHeavy[parentOriginal], let child = originalToHeavy[childOriginal] else { return nil }
        var mask: UInt64 = 0
        for serial in branch.descendants {
            if let original = serialToOriginal[serial], let heavy = originalToHeavy[original] { mask |= UInt64(1) << UInt64(heavy) }
        }
        return TorsionData(parent: UInt32(parent), child: UInt32(child), maskLow: UInt32(mask & 0xffffffff), maskHigh: UInt32(mask >> 32))
    }
}

func buildPairs(typed: TypedLigand, torsions: [TorsionData]) -> [PairData] {
    func graphDistance(_ start: Int, _ goal: Int, limit: Int = 3) -> Int? {
        var visited: Set<Int> = [start]
        var frontier = [start]
        for distance in 1...limit {
            frontier = frontier.flatMap { typed.bonds[$0] }.filter { visited.insert($0).inserted }
            if frontier.contains(goal) { return distance }
        }
        return nil
    }
    func inMask(_ heavy: Int, _ torsion: TorsionData) -> Bool {
        heavy < 32 ? (torsion.maskLow & (UInt32(1) << UInt32(heavy))) != 0 : (torsion.maskHigh & (UInt32(1) << UInt32(heavy - 32))) != 0
    }
    var pairs: [PairData] = []
    for a in typed.atoms.indices {
        for b in typed.atoms.indices where b > a {
            let originalA = typed.originalIndices[a], originalB = typed.originalIndices[b]
            if graphDistance(originalA, originalB) != nil { continue }
            if !torsions.contains(where: { inMask(a, $0) != inMask(b, $0) }) { continue }
            pairs.append(PairData(a: UInt32(a), b: UInt32(b), typeA: typed.xsTypes[a], typeB: typed.xsTypes[b]))
        }
    }
    return pairs
}

func extractExactPairs(extractor: URL, ligandURL: URL, typed: TypedLigand) throws -> [PairData] {
    let process = Process()
    let pipe = Pipe()
    process.executableURL = extractor
    process.arguments = [ligandURL.path]
    process.standardOutput = pipe
    process.standardError = FileHandle.standardError
    try process.run()
    let data = pipe.fileHandleForReading.readDataToEndOfFile()
    process.waitUntilExit()
    guard process.terminationStatus == 0 else { throw NSError(domain: "VinaMetal", code: 4) }
    let lines = String(decoding: data, as: UTF8.self).split(whereSeparator: \.isNewline)
    var vinaToHeavy: [Int: Int] = [:]
    var vinaTypes: [Int: UInt32] = [:]
    for line in lines where line.hasPrefix("ATOM ") {
        let fields = line.split(whereSeparator: \.isWhitespace)
        guard fields.count == 6, let vinaIndex = Int(fields[1]), let xs = UInt32(fields[2]),
              let x = Float(fields[3]), let y = Float(fields[4]), let z = Float(fields[5]), xs < 32 else { continue }
        let best = typed.atoms.indices.min { lhs, rhs in
            let a=typed.atoms[lhs], b=typed.atoms[rhs]
            let da=(a.x-x)*(a.x-x)+(a.y-y)*(a.y-y)+(a.z-z)*(a.z-z)
            let db=(b.x-x)*(b.x-x)+(b.y-y)*(b.y-y)+(b.z-z)*(b.z-z)
            return da < db
        }!
        vinaToHeavy[vinaIndex] = best
        vinaTypes[vinaIndex] = xs
    }
    var pairs: [PairData] = []
    for line in lines where line.hasPrefix("PAIR ") {
        let fields = line.split(whereSeparator: \.isWhitespace)
        guard fields.count == 3, let vinaA = Int(fields[1]), let vinaB = Int(fields[2]),
              let a = vinaToHeavy[vinaA], let b = vinaToHeavy[vinaB],
              let typeA = vinaTypes[vinaA], let typeB = vinaTypes[vinaB] else { continue }
        pairs.append(PairData(a: UInt32(a), b: UInt32(b), typeA: typeA, typeB: typeB))
    }
    return pairs
}

func makeTranslations(count: Int) -> [PackedFloat4] {
    var state: UInt64 = 0x9E3779B97F4A7C15
    func random() -> Float {
        state = state &* 6364136223846793005 &+ 1442695040888963407
        return Float((state >> 40) & 0xFFFFFF) / Float(0xFFFFFF)
    }
    return (0..<count).map { _ in
        PackedFloat4(x: (random() - 0.5) * 2, y: (random() - 0.5) * 2, z: (random() - 0.5) * 2, w: 0)
    }
}

@inline(__always)
func cpuScore(maps: [Float], atoms: [PackedFloat4], translation: PackedFloat4, params: GridParams) -> PackedFloat4 {
    let dx = Int(params.dimX), dy = Int(params.dimY), dz = Int(params.dimZ)
    let stride = dx * dy * dz
    var energy: Float = 0
    var totalGradient: [Float] = [0, 0, 0]
    for atom in atoms {
        var s = [
            (atom.x + translation.x - params.originX) / params.spacing,
            (atom.y + translation.y - params.originY) / params.spacing,
            (atom.z + translation.z - params.originZ) / params.spacing,
        ]
        var miss: [Float] = [0, 0, 0]
        var region: [Float] = [0, 0, 0]
        let dims = [dx, dy, dz]
        for axis in 0..<3 {
            if s[axis] < 0 { miss[axis] = -s[axis]; region[axis] = -1; s[axis] = 0 }
            else if s[axis] >= Float(dims[axis] - 1) { miss[axis] = s[axis] - Float(dims[axis] - 1); region[axis] = 1; s[axis] = Float(dims[axis] - 1) }
        }
        let x0 = min(Int(s[0]), dx - 2), y0 = min(Int(s[1]), dy - 2), z0 = min(Int(s[2]), dz - 2)
        let x = min(max(s[0] - Float(x0), 0), 1), y = min(max(s[1] - Float(y0), 0), 1), z = min(max(s[2] - Float(z0), 0), 1)
        let mx = 1 - x, my = 1 - y, mz = 1 - z
        let base = Int(atom.w) * stride
        func value(_ xi: Int, _ yi: Int, _ zi: Int) -> Float { maps[base + xi + dx * (yi + dy * zi)] }
        let f000 = value(x0, y0, z0), f100 = value(x0 + 1, y0, z0)
        let f010 = value(x0, y0 + 1, z0), f110 = value(x0 + 1, y0 + 1, z0)
        let f001 = value(x0, y0, z0 + 1), f101 = value(x0 + 1, y0, z0 + 1)
        let f011 = value(x0, y0 + 1, z0 + 1), f111 = value(x0 + 1, y0 + 1, z0 + 1)
        var interpolated = f000 * mx * my * mz + f100 * x * my * mz
            + value(x0, y0 + 1, z0) * mx * y * mz + value(x0 + 1, y0 + 1, z0) * x * y * mz
            + value(x0, y0, z0 + 1) * mx * my * z + value(x0 + 1, y0, z0 + 1) * x * my * z
            + value(x0, y0 + 1, z0 + 1) * mx * y * z + value(x0 + 1, y0 + 1, z0 + 1) * x * y * z
        var gradient = [
            f000 * -my * mz + f100 * my * mz + f010 * -y * mz + f110 * y * mz + f001 * -my * z + f101 * my * z + f011 * -y * z + f111 * y * z,
            f000 * mx * -mz + f100 * x * -mz + f010 * mx * mz + f110 * x * mz + f001 * mx * -z + f101 * x * -z + f011 * mx * z + f111 * x * z,
            f000 * mx * my * -1 + f100 * x * my * -1 + f010 * mx * y * -1 + f110 * x * y * -1 + f001 * mx * my + f101 * x * my + f011 * mx * y + f111 * x * y,
        ]
        if interpolated > 0 {
            let curl = params.curlV / (params.curlV + interpolated)
            interpolated *= curl
            gradient = gradient.map { $0 * curl * curl }
        }
        energy += interpolated + params.slope * (miss[0] + miss[1] + miss[2]) * params.spacing
        for axis in 0..<3 { totalGradient[axis] += gradient[axis] / params.spacing + params.slope * region[axis] }
    }
    return PackedFloat4(x: totalGradient[0], y: totalGradient[1], z: totalGradient[2], w: energy)
}

func cpuIntraEnergy(atoms: [PackedFloat4], pairs: [PairData]) -> Float {
    let radii: [Float] = [1.9,1.9,1.8,1.8,1.8,1.8,1.7,1.7,1.7,1.7,2.0,2.1,1.5,1.8,2.0,2.2]
    func hydrophobic(_ type: UInt32) -> Bool { [0,12,13,14,15].contains(type) }
    func acceptor(_ type: UInt32) -> Bool { [4,5,8,9].contains(type) }
    func donor(_ type: UInt32) -> Bool { [3,5,7,9,18].contains(type) }
    func slope(_ bad: Float, _ good: Float, _ value: Float) -> Float {
        if bad < good { if value <= bad { return 0 }; if value >= good { return 1 } }
        else { if value >= bad { return 0 }; if value <= good { return 1 } }
        return (value - bad) / (good - bad)
    }
    var energy: Float = 0
    for pair in pairs {
        let a = atoms[Int(pair.a)], b = atoms[Int(pair.b)]
        let dx = a.x-b.x, dy = a.y-b.y, dz = a.z-b.z
        let distance = sqrt(dx*dx+dy*dy+dz*dz)
        if distance >= 8 { continue }
        let d = distance-radii[min(Int(pair.typeA),15)]-radii[min(Int(pair.typeB),15)]
        let gauss1=exp(-pow(d/0.5,2)), gauss2=exp(-pow((d-3)/2,2)), repulsion=d<0 ? d*d : 0
        let hydro=hydrophobic(pair.typeA)&&hydrophobic(pair.typeB) ? slope(1.5,0.5,d) : 0
        let hbond=(donor(pair.typeA)&&acceptor(pair.typeB))||(donor(pair.typeB)&&acceptor(pair.typeA)) ? slope(0,-0.7,d) : 0
        energy += -0.035579*gauss1-0.005156*gauss2+0.840245*repulsion-0.035069*hydro-0.587439*hbond
    }
    return energy
}

let root = URL(fileURLWithPath: FileManager.default.currentDirectoryPath)
let mapsDirectory = URL(fileURLWithPath: stringArgument("--maps", default: root.appending(path: "maps").path))
let mapURLs = try FileManager.default.contentsOfDirectory(at: mapsDirectory, includingPropertiesForKeys: nil)
    .filter { $0.pathExtension == "map" }.sorted { $0.lastPathComponent < $1.lastPathComponent }
guard !mapURLs.isEmpty else { fatalError("No .map files in \(mapsDirectory.path)") }
let parsedMaps = try mapURLs.map(parseMap)
let mapNames = mapURLs.map { $0.deletingPathExtension().pathExtension }
let mapIndices = Dictionary(uniqueKeysWithValues: mapNames.enumerated().map { ($0.element, $0.offset) })
let first = parsedMaps[0]
guard parsedMaps.allSatisfy({ $0.dims == first.dims && $0.spacing == first.spacing }) else { fatalError("Inconsistent map geometry") }
let maps = parsedMaps.flatMap(\.values)
let ligandURL = URL(fileURLWithPath: stringArgument("--ligand", default: root.deletingLastPathComponent().appending(path: "data/1iep_ligand.pdbqt").path))
let ligand = try parseLigand(ligandURL)
let typedLigand = try typedHeavyAtoms(ligand.atoms, mapIndices: mapIndices)
let atoms = typedLigand.atoms
let torsions = parseTorsions(lines: ligand.lines, ligand: ligand.atoms, heavyOriginalIndices: typedLigand.originalIndices)
let inferredPairs = buildPairs(typed: typedLigand, torsions: torsions)
let extractorURL = root.appending(path: "tools/extract_vina_model")
let interactingPairs = FileManager.default.isExecutableFile(atPath: extractorURL.path)
    ? try extractExactPairs(extractor: extractorURL, ligandURL: ligandURL, typed: typedLigand)
    : inferredPairs
let origin = (
    first.center.0 - Float(first.dims.0 - 1) * first.spacing / 2,
    first.center.1 - Float(first.dims.1 - 1) * first.spacing / 2,
    first.center.2 - Float(first.dims.2 - 1) * first.spacing / 2
)

guard let device = MTLCreateSystemDefaultDevice(), let queue = device.makeCommandQueue() else { fatalError("No Metal device") }
let shaderURL = Bundle.module.url(forResource: "grid_score", withExtension: "metal", subdirectory: "Shaders")!
let shader = try String(contentsOf: shaderURL, encoding: .utf8)
let library = try device.makeLibrary(source: shader, options: nil)
let pipeline = try device.makeComputePipelineState(function: library.makeFunction(name: "scorePoses")!)

func buffer<T>(_ values: [T], options: MTLResourceOptions = .storageModeShared) -> MTLBuffer {
    values.withUnsafeBytes { bytes in
        device.makeBuffer(bytes: bytes.baseAddress!, length: bytes.count, options: options)!
    }
}

func quaternionRotate(_ q: PackedFloat4, _ v: PackedFloat4) -> PackedFloat4 {
    let qx = q.x, qy = q.y, qz = q.z, qw = q.w
    let tx = 2 * (qy * v.z - qz * v.y)
    let ty = 2 * (qz * v.x - qx * v.z)
    let tz = 2 * (qx * v.y - qy * v.x)
    return PackedFloat4(
        x: v.x + qw * tx + (qy * tz - qz * ty),
        y: v.y + qw * ty + (qz * tx - qx * tz),
        z: v.z + qw * tz + (qx * ty - qy * tx),
        w: v.w
    )
}

func replaceCoordinates(in line: String, with position: PackedFloat4) -> String {
    var bytes = Array(line.utf8)
    if bytes.count < 54 { bytes += Array(repeating: 32, count: 54 - bytes.count) }
    let coordinates = String(format: "%8.3f%8.3f%8.3f", position.x, position.y, position.z)
    bytes.replaceSubrange(30..<54, with: coordinates.utf8)
    return String(decoding: bytes, as: UTF8.self)
}

func applyTorsionsToAll(lines: [String], atoms: [LigandAtom], centered: [PackedFloat4], angles: [Float]) -> [PackedFloat4] {
    struct Branch { var parent: Int; var child: Int; var descendants: Set<Int> }
    var branches: [Branch] = [], stack: [Int] = []
    for line in lines {
        let fields = line.split(whereSeparator: \.isWhitespace)
        if fields.first == "BRANCH", fields.count >= 3, let parent = Int(fields[1]), let child = Int(fields[2]) {
            branches.append(Branch(parent: parent, child: child, descendants: [])); stack.append(branches.count - 1)
        } else if fields.first == "ENDBRANCH" { _ = stack.popLast() }
        else if (fields.first == "ATOM" || fields.first == "HETATM"), fields.count > 1, let serial = Int(fields[1]) {
            for index in stack { branches[index].descendants.insert(serial) }
        }
    }
    let serialToIndex = Dictionary(uniqueKeysWithValues: atoms.indices.map { (atoms[$0].serial, $0) })
    var coordinates = centered
    for (torsionIndex, branch) in branches.enumerated() where torsionIndex < angles.count {
        guard let parent = serialToIndex[branch.parent], let child = serialToIndex[branch.child] else { continue }
        let axisPoint = coordinates[parent]
        var ax = coordinates[child].x - axisPoint.x, ay = coordinates[child].y - axisPoint.y, az = coordinates[child].z - axisPoint.z
        let length = sqrt(ax * ax + ay * ay + az * az); ax /= length; ay /= length; az /= length
        let c = cos(angles[torsionIndex]), s = sin(angles[torsionIndex])
        for serial in branch.descendants {
            guard let index = serialToIndex[serial] else { continue }
            let rx = coordinates[index].x - axisPoint.x, ry = coordinates[index].y - axisPoint.y, rz = coordinates[index].z - axisPoint.z
            let crossX = ay * rz - az * ry, crossY = az * rx - ax * rz, crossZ = ax * ry - ay * rx
            let dot = ax * rx + ay * ry + az * rz
            coordinates[index] = PackedFloat4(
                x: axisPoint.x + rx * c + crossX * s + ax * dot * (1 - c),
                y: axisPoint.y + ry * c + crossY * s + ay * dot * (1 - c),
                z: axisPoint.z + rz * c + crossZ * s + az * dot * (1 - c), w: 0)
        }
    }
    return coordinates
}

func integerArgument(_ name: String, default defaultValue: Int) -> Int {
    guard let index = CommandLine.arguments.firstIndex(of: name), index + 1 < CommandLine.arguments.count else { return defaultValue }
    return Int(CommandLine.arguments[index + 1]) ?? defaultValue
}

func stringArgument(_ name: String, default defaultValue: String) -> String {
    guard let index = CommandLine.arguments.firstIndex(of: name), index + 1 < CommandLine.arguments.count else { return defaultValue }
    return CommandLine.arguments[index + 1]
}

func floatArgument(_ name: String, default defaultValue: Float) -> Float {
    guard let index = CommandLine.arguments.firstIndex(of: name), index + 1 < CommandLine.arguments.count else { return defaultValue }
    return Float(CommandLine.arguments[index + 1]) ?? defaultValue
}

func hasArgument(_ name: String) -> Bool { CommandLine.arguments.contains(name) }

func milliseconds(_ duration: Duration) -> Double {
    let components = duration.components
    return Double(components.seconds) * 1_000 + Double(components.attoseconds) / 1e15
}

let mapsBuffer = buffer(maps)
let atomsBuffer = buffer(atoms)
print("Device: \(device.name), unified memory: \(device.hasUnifiedMemory)")
print("Maps: \(parsedMaps.count) × \(first.dims.0)×\(first.dims.1)×\(first.dims.2); atoms: \(atoms.count)")
let crystalParams = GridParams(dimX: UInt32(first.dims.0), dimY: UInt32(first.dims.1), dimZ: UInt32(first.dims.2), mapCount: UInt32(parsedMaps.count), atomCount: UInt32(atoms.count), poseCount: 1, spacing: first.spacing, slope: 1e6, originX: origin.0, originY: origin.1, originZ: origin.2, curlV: 1000)
let crystalGridEnergy = cpuScore(maps: maps, atoms: atoms, translation: PackedFloat4(x: 0, y: 0, z: 0, w: 0), params: crystalParams).w
print("Input-pose grid energy: \(String(format: "%.3f", crystalGridEnergy)) kcal/mol")
print("Input-pose internal energy: \(String(format: "%.3f", cpuIntraEnergy(atoms: atoms, pairs: interactingPairs))) kcal/mol")
if CommandLine.arguments.contains("--dock") {
    let flexible = CommandLine.arguments.contains("--flexible")
    let exhaustiveness = integerArgument("--exhaustiveness", default: 8)
    let degreesOfFreedom = 6 + torsions.count
    let heuristic = ligand.atoms.count + 10 * degreesOfFreedom
    let vinaGlobalSteps = 105 * (50 + heuristic)
    let vinaLocalSteps = (25 + ligand.atoms.count) / 3
    let totalMutations = exhaustiveness * vinaGlobalSteps
    let laneTarget = min(8_192, max(256, totalMutations / 32))
    var adaptiveLanes = 256
    while adaptiveLanes * 2 <= laneTarget { adaptiveLanes *= 2 }
    let lanes = integerArgument("--lanes", default: adaptiveLanes)
    let steps = integerArgument("--steps", default: Int(ceil(Double(totalMutations) / Double(lanes))))
    let localSteps = integerArgument("--local-steps", default: vinaLocalSteps)
    let center = ligand.atoms.reduce(PackedFloat4(x: 0, y: 0, z: 0, w: 0)) { partial, atom in
        PackedFloat4(x: partial.x + atom.position.x, y: partial.y + atom.position.y, z: partial.z + atom.position.z, w: 0)
    }
    let inverseCount = 1 / Float(ligand.atoms.count)
    let ligandCenter = PackedFloat4(x: center.x * inverseCount, y: center.y * inverseCount, z: center.z * inverseCount, w: 0)
    let localAtoms = atoms.map { atom in
        PackedFloat4(x: atom.x - ligandCenter.x, y: atom.y - ligandCenter.y, z: atom.z - ligandCenter.z, w: atom.w)
    }
    let gyrationRadius = sqrt(localAtoms.reduce(Float(0)) { $0 + $1.x*$1.x + $1.y*$1.y + $1.z*$1.z } / Float(max(localAtoms.count, 1)))
    let localAtomsBuffer = buffer(localAtoms)
    let poseBuffer = device.makeBuffer(length: lanes * MemoryLayout<PackedFloat4>.stride, options: .storageModeShared)!
    let orientationBuffer = device.makeBuffer(length: lanes * MemoryLayout<PackedFloat4>.stride, options: .storageModeShared)!
    let angleBuffer = device.makeBuffer(length: lanes * 2 * MemoryLayout<PackedFloat4>.stride, options: .storageModeShared)!
    let torsionBuffer = buffer(torsions)
    let pairBuffer = buffer(interactingPairs)
    var grid = GridParams(dimX: UInt32(first.dims.0), dimY: UInt32(first.dims.1), dimZ: UInt32(first.dims.2), mapCount: UInt32(parsedMaps.count), atomCount: UInt32(localAtoms.count), poseCount: UInt32(lanes), spacing: first.spacing, slope: 1e6, originX: origin.0, originY: origin.1, originZ: origin.2, curlV: 10)
    var search = SearchParams(
        lanes: UInt32(lanes), steps: UInt32(steps), seed: UInt32(integerArgument("--seed", default: 20_260_717)), localSteps: UInt32(localSteps),
        torsionCount: UInt32(torsions.count), pairCount: UInt32(interactingPairs.count),
        centerX: first.center.0, centerY: first.center.1, centerZ: first.center.2,
        spanX: floatArgument("--span-x", default: Float(first.dims.0 - 1) * first.spacing),
        spanY: floatArgument("--span-y", default: Float(first.dims.1 - 1) * first.spacing),
        spanZ: floatArgument("--span-z", default: Float(first.dims.2 - 1) * first.spacing), translationMutation: 2.0,
        rotationMutation: gyrationRadius > 1e-6 ? 2.0 / gyrationRadius : 0, temperature: 1.2, gradientStep: 0.3
    )
    let dockPipeline = try device.makeComputePipelineState(function: library.makeFunction(name: flexible ? "flexibleDock" : "rigidDock")!)
    let started = ContinuousClock.now
    let commandBuffer = queue.makeCommandBuffer()!
    let encoder = commandBuffer.makeComputeCommandEncoder()!
    encoder.setComputePipelineState(dockPipeline)
    encoder.setBuffer(mapsBuffer, offset: 0, index: 0)
    encoder.setBuffer(localAtomsBuffer, offset: 0, index: 1)
    encoder.setBuffer(poseBuffer, offset: 0, index: 2)
    encoder.setBuffer(orientationBuffer, offset: 0, index: 3)
    encoder.setBytes(&grid, length: MemoryLayout<GridParams>.stride, index: 4)
    encoder.setBytes(&search, length: MemoryLayout<SearchParams>.stride, index: 5)
    if flexible {
        encoder.setBuffer(torsionBuffer, offset: 0, index: 6)
        encoder.setBuffer(pairBuffer, offset: 0, index: 7)
        encoder.setBuffer(angleBuffer, offset: 0, index: 8)
    }
    let width = min(dockPipeline.threadExecutionWidth * 4, dockPipeline.maxTotalThreadsPerThreadgroup)
    encoder.dispatchThreads(MTLSize(width: lanes, height: 1, depth: 1), threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1))
    encoder.endEncoding()
    commandBuffer.commit()
    commandBuffer.waitUntilCompleted()
    if let error = commandBuffer.error { fatalError("Metal docking failed: \(error)") }
    let elapsedMs = milliseconds(started.duration(to: .now))
    let poses = UnsafeBufferPointer(start: poseBuffer.contents().bindMemory(to: PackedFloat4.self, capacity: lanes), count: lanes)
    let orientations = UnsafeBufferPointer(start: orientationBuffer.contents().bindMemory(to: PackedFloat4.self, capacity: lanes), count: lanes)
    let angleResults = UnsafeBufferPointer(start: angleBuffer.contents().bindMemory(to: PackedFloat4.self, capacity: lanes * 2), count: lanes * 2)
    let centeredAllAtoms = ligand.atoms.map { PackedFloat4(x: $0.position.x - ligandCenter.x, y: $0.position.y - ligandCenter.y, z: $0.position.z - ligandCenter.z, w: 0) }
    func laneAngles(_ index: Int) -> [Float] {
        guard flexible else { return Array(repeating: 0, count: 8) }
        let a = angleResults[index * 2], b = angleResults[index * 2 + 1]
        return [a.x, a.y, a.z, a.w, b.x, b.y, b.z, b.w]
    }
    func laneHeavyCoordinates(_ index: Int) -> [PackedFloat4] {
        let torsioned = flexible ? applyTorsionsToAll(lines: ligand.lines, atoms: ligand.atoms, centered: centeredAllAtoms, angles: laneAngles(index)) : centeredAllAtoms
        return ligand.atoms.indices.compactMap { atomIndex in
            guard ligand.atoms[atomIndex].adType != "H" && ligand.atoms[atomIndex].adType != "HD" else { return nil }
            let rotated = quaternionRotate(orientations[index], torsioned[atomIndex])
            return PackedFloat4(x: rotated.x + poses[index].x, y: rotated.y + poses[index].y, z: rotated.z + poses[index].z, w: 0)
        }
    }
    func coordinateRMSD(_ lhs: [PackedFloat4], _ rhs: [PackedFloat4]) -> Float {
        guard lhs.count == rhs.count && !lhs.isEmpty else { return .infinity }
        let squared = zip(lhs, rhs).reduce(Float(0)) { sum, pair in
            let dx=pair.0.x-pair.1.x, dy=pair.0.y-pair.1.y, dz=pair.0.z-pair.1.z
            return sum + dx*dx + dy*dy + dz*dz
        }
        return sqrt(squared / Float(lhs.count))
    }
    let numModes = integerArgument("--num-modes", default: 9)
    let minRMSD = floatArgument("--min-rmsd", default: 1.0)
    var clustered: [(index: Int, coordinates: [PackedFloat4])] = []
    for index in poses.indices.sorted(by: { poses[$0].w < poses[$1].w }) {
        let coordinates = laneHeavyCoordinates(index)
        if clustered.allSatisfy({ coordinateRMSD(coordinates, $0.coordinates) >= minRMSD }) {
            clustered.append((index, coordinates))
            if clustered.count == numModes { break }
        }
    }
    let bestIndex = clustered.first!.index
    let bestPose = poses[bestIndex]
    let bestOrientation = orientations[bestIndex]
    let bestAngles = laneAngles(bestIndex)
    let torsionedAllAtoms = flexible ? applyTorsionsToAll(lines: ligand.lines, atoms: ligand.atoms, centered: centeredAllAtoms, angles: bestAngles) : centeredAllAtoms
    var outputLines = ligand.lines
    for (atomIndex, atom) in ligand.atoms.enumerated() {
        let rotated = quaternionRotate(bestOrientation, torsionedAllAtoms[atomIndex])
        let transformed = PackedFloat4(x: rotated.x + bestPose.x, y: rotated.y + bestPose.y, z: rotated.z + bestPose.z, w: 0)
        outputLines[atom.lineIndex] = replaceCoordinates(in: outputLines[atom.lineIndex], with: transformed)
    }
    let resultsDirectory = root.deletingLastPathComponent().appending(path: "results/metal")
    try FileManager.default.createDirectory(at: resultsDirectory, withIntermediateDirectories: true)
    let defaultOutput = resultsDirectory.appending(path: flexible ? "1iep_metal_flexible_out.pdbqt" : "1iep_metal_rigid_out.pdbqt").path
    let outputURL = URL(fileURLWithPath: stringArgument("--output", default: defaultOutput))
    try FileManager.default.createDirectory(at: outputURL.deletingLastPathComponent(), withIntermediateDirectories: true)
    try (outputLines.joined(separator: "\n") + "\n").write(to: outputURL, atomically: true, encoding: .utf8)
    let modesURL = URL(fileURLWithPath: stringArgument("--modes-output", default: outputURL.deletingPathExtension().path + "_modes.pdbqt"))
    var modeBlocks: [String] = []
    var clusteredRawURLs: [URL] = []
    for (mode, candidate) in clustered.enumerated() {
        let index = candidate.index, pose = poses[index], orientation = orientations[index]
        let torsioned = flexible ? applyTorsionsToAll(lines: ligand.lines, atoms: ligand.atoms, centered: centeredAllAtoms, angles: laneAngles(index)) : centeredAllAtoms
        var lines = ligand.lines
        for (atomIndex, atom) in ligand.atoms.enumerated() {
            let rotated = quaternionRotate(orientation, torsioned[atomIndex])
            lines[atom.lineIndex] = replaceCoordinates(in: lines[atom.lineIndex], with: PackedFloat4(x: rotated.x+pose.x,y: rotated.y+pose.y,z: rotated.z+pose.z,w: 0))
        }
        modeBlocks.append("MODEL \(mode + 1)\nREMARK METAL SEARCH ENERGY: \(String(format: "%.6f", pose.w))\n" + lines.joined(separator: "\n") + "\nENDMDL")
        let candidateURL = URL(fileURLWithPath: outputURL.deletingPathExtension().path + "_mode\(mode + 1).pdbqt")
        try (lines.joined(separator: "\n") + "\n").write(to: candidateURL, atomically: true, encoding: .utf8)
        clusteredRawURLs.append(candidateURL)
    }
    try (modeBlocks.joined(separator: "\n") + "\n").write(to: modesURL, atomically: true, encoding: .utf8)
    print("Metal \(flexible ? "flexible" : "rigid") docking: \(lanes) lanes × \(steps) steps, \(String(format: "%.1f", elapsedMs)) ms")
    print("Vina-compatible effort: exhaustiveness \(exhaustiveness), heuristic global steps \(vinaGlobalSteps), local steps \(localSteps)")
    print("Torsions: \(torsions.count); intramolecular pairs: \(interactingPairs.count)")
    print("Clustered minima: \(clustered.count) at \(String(format: "%.2f", minRMSD)) Å RMSD")
    print("Best Metal search energy: \(String(format: "%.3f", bestPose.w)) kcal/mol")
    print("Output: \(outputURL.path)")
    print("Clustered modes output: \(modesURL.path)")
    if hasArgument("--vina-receptor") && hasArgument("--vina-config") {
        let vinaBinary = URL(fileURLWithPath: stringArgument("--vina-binary", default: root.deletingLastPathComponent().appending(path: "bin/vina").path))
        let receptor = stringArgument("--vina-receptor", default: "")
        _ = stringArgument("--vina-config", default: "")
        let finalizedURL = URL(fileURLWithPath: stringArgument("--vina-output", default: outputURL.deletingPathExtension().path + "_vina_final.pdbqt"))
        let finalizationStart = ContinuousClock.now
        var jobs: [(process: Process, pipe: Pipe, output: URL)] = []
        for (mode, candidateURL) in clusteredRawURLs.enumerated() {
            let candidateOutput = URL(fileURLWithPath: finalizedURL.deletingPathExtension().path + "_mode\(mode + 1).pdbqt")
            let process = Process(), pipe = Pipe()
            process.executableURL = vinaBinary
            // The Metal grid already enforced the requested docking box. Autoboxing
            // the local refinement prevents Vina from rejecting a valid candidate
            // merely because an outer ligand atom lies on that box's boundary.
            process.arguments = ["--receptor", receptor, "--ligand", candidateURL.path, "--autobox", "--local_only", "--out", candidateOutput.path]
            process.standardOutput = pipe; process.standardError = pipe
            try process.run()
            jobs.append((process, pipe, candidateOutput))
        }
        let scorePattern = try NSRegularExpression(pattern: #"Estimated Free Energy of Binding\s*:\s*(-?[0-9.]+)"#)
        var finalized: [(score: Double, output: URL)] = []
        for job in jobs {
            let data = job.pipe.fileHandleForReading.readDataToEndOfFile()
            job.process.waitUntilExit()
            guard job.process.terminationStatus == 0 else {
                throw NSError(domain: "VinaMetal", code: 5, userInfo: [NSLocalizedDescriptionKey: String(decoding: data, as: UTF8.self)])
            }
            let log = String(decoding: data, as: UTF8.self), range = NSRange(log.startIndex..<log.endIndex, in: log)
            guard let match = scorePattern.firstMatch(in: log, range: range), let scoreRange = Range(match.range(at: 1), in: log), let score = Double(log[scoreRange]) else {
                throw NSError(domain: "VinaMetal", code: 6, userInfo: [NSLocalizedDescriptionKey: "Could not parse authoritative Vina score"])
            }
            finalized.append((score, job.output))
        }
        let bestFinal = finalized.min(by: { $0.score < $1.score })!
        try Data(contentsOf: bestFinal.output).write(to: finalizedURL, options: .atomic)
        let finalizationMs = milliseconds(finalizationStart.duration(to: .now))
        print("Authoritative Vina 1.2.7 score: \(String(format: "%.3f", bestFinal.score)) kcal/mol")
        print("Authoritative modes refined: \(finalized.count)")
        print("Vina finalization: \(String(format: "%.1f", finalizationMs)) ms")
        print("Finalized output: \(finalizedURL.path)")
    }
    exit(0)
}

print("poses,cpu_ms,gpu_ms,speedup,max_energy_error,p99_gradient_error,max_gradient_error")
for poseCount in [1, 64, 1_024, 16_384, 65_536] {
    let translations = makeTranslations(count: poseCount)
    let translationsBuffer = buffer(translations)
    let outputBuffer = device.makeBuffer(length: poseCount * MemoryLayout<PackedFloat4>.stride, options: .storageModeShared)!
    var params = GridParams(dimX: UInt32(first.dims.0), dimY: UInt32(first.dims.1), dimZ: UInt32(first.dims.2), mapCount: UInt32(parsedMaps.count), atomCount: UInt32(atoms.count), poseCount: UInt32(poseCount), spacing: first.spacing, slope: 1e6, originX: origin.0, originY: origin.1, originZ: origin.2, curlV: 10)

    let cpuStart = ContinuousClock.now
    let cpu = translations.map { cpuScore(maps: maps, atoms: atoms, translation: $0, params: params) }
    let cpuMs = milliseconds(cpuStart.duration(to: .now))

    func runGPU() {
        let commandBuffer = queue.makeCommandBuffer()!
        let encoder = commandBuffer.makeComputeCommandEncoder()!
        encoder.setComputePipelineState(pipeline)
        encoder.setBuffer(mapsBuffer, offset: 0, index: 0)
        encoder.setBuffer(atomsBuffer, offset: 0, index: 1)
        encoder.setBuffer(translationsBuffer, offset: 0, index: 2)
        encoder.setBuffer(outputBuffer, offset: 0, index: 3)
        encoder.setBytes(&params, length: MemoryLayout<GridParams>.stride, index: 4)
        let width = min(pipeline.threadExecutionWidth * 4, pipeline.maxTotalThreadsPerThreadgroup)
        encoder.dispatchThreads(MTLSize(width: poseCount, height: 1, depth: 1), threadsPerThreadgroup: MTLSize(width: width, height: 1, depth: 1))
        encoder.endEncoding()
        commandBuffer.commit()
        commandBuffer.waitUntilCompleted()
        if let error = commandBuffer.error { fatalError("Metal command failed: \(error)") }
    }

    runGPU()
    let gpuStart = ContinuousClock.now
    runGPU()
    let gpuMs = milliseconds(gpuStart.duration(to: .now))
    let gpu = Array(UnsafeBufferPointer(start: outputBuffer.contents().bindMemory(to: PackedFloat4.self, capacity: poseCount), count: poseCount))
    let maxEnergyError = zip(cpu, gpu).map { abs($0.w - $1.w) }.max() ?? 0
    let gradientErrors = zip(cpu, gpu).map { max(abs($0.x - $1.x), abs($0.y - $1.y), abs($0.z - $1.z)) }.sorted()
    let maxGradientError = gradientErrors.last ?? 0
    let p99GradientError = gradientErrors[min(Int(Double(gradientErrors.count) * 0.99), gradientErrors.count - 1)]
    guard maxEnergyError <= 0.01 && p99GradientError <= 0.01 else {
        fatalError("CPU/Metal mismatch: energy \(maxEnergyError), p99 gradient \(p99GradientError), max gradient \(maxGradientError)")
    }
    print("\(poseCount),\(String(format: "%.3f", cpuMs)),\(String(format: "%.3f", gpuMs)),\(String(format: "%.2f", cpuMs / gpuMs)),\(String(format: "%.6f", maxEnergyError)),\(String(format: "%.6f", p99GradientError)),\(String(format: "%.6f", maxGradientError))")
}
