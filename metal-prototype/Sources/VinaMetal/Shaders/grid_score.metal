#include <metal_stdlib>
using namespace metal;

struct GridParams {
    uint dimX;
    uint dimY;
    uint dimZ;
    uint mapCount;
    uint atomCount;
    uint poseCount;
    float spacing;
    float slope;
    float originX;
    float originY;
    float originZ;
    float curlV;
};

inline uint gridIndex(uint x, uint y, uint z, constant GridParams& p) {
    return x + p.dimX * (y + p.dimY * z);
}

kernel void scorePoses(
    device const float* maps [[buffer(0)]],
    device const packed_float4* atoms [[buffer(1)]],
    device const packed_float4* translations [[buffer(2)]],
    device packed_float4* scores [[buffer(3)]],
    constant GridParams& params [[buffer(4)]],
    uint pose [[thread_position_in_grid]]
) {
    if (pose >= params.poseCount) return;

    const float3 translation = translations[pose].xyz;
    const uint mapStride = params.dimX * params.dimY * params.dimZ;
    float energy = 0.0f;
    float3 totalGradient = 0.0f;

    for (uint atomIndex = 0; atomIndex < params.atomCount; ++atomIndex) {
        const packed_float4 atom = atoms[atomIndex];
        const float3 location = atom.xyz + translation;
        const uint mapIndex = min(uint(atom.w), params.mapCount - 1);
        float3 s = (location - float3(params.originX, params.originY, params.originZ)) / params.spacing;
        float3 miss = 0.0f;
        int3 region = 0;

        if (s.x < 0.0f) { miss.x = -s.x; region.x = -1; s.x = 0.0f; }
        else if (s.x >= float(params.dimX - 1)) { miss.x = s.x - float(params.dimX - 1); region.x = 1; s.x = float(params.dimX - 1); }
        if (s.y < 0.0f) { miss.y = -s.y; region.y = -1; s.y = 0.0f; }
        else if (s.y >= float(params.dimY - 1)) { miss.y = s.y - float(params.dimY - 1); region.y = 1; s.y = float(params.dimY - 1); }
        if (s.z < 0.0f) { miss.z = -s.z; region.z = -1; s.z = 0.0f; }
        else if (s.z >= float(params.dimZ - 1)) { miss.z = s.z - float(params.dimZ - 1); region.z = 1; s.z = float(params.dimZ - 1); }

        const uint x0 = min(uint(s.x), params.dimX - 2);
        const uint y0 = min(uint(s.y), params.dimY - 2);
        const uint z0 = min(uint(s.z), params.dimZ - 2);
        const float3 fraction = clamp(s - float3(x0, y0, z0), 0.0f, 1.0f);
        const float3 inverse = 1.0f - fraction;
        const uint base = mapIndex * mapStride;

        const float f000 = maps[base + gridIndex(x0,     y0,     z0,     params)];
        const float f100 = maps[base + gridIndex(x0 + 1, y0,     z0,     params)];
        const float f010 = maps[base + gridIndex(x0,     y0 + 1, z0,     params)];
        const float f110 = maps[base + gridIndex(x0 + 1, y0 + 1, z0,     params)];
        const float f001 = maps[base + gridIndex(x0,     y0,     z0 + 1, params)];
        const float f101 = maps[base + gridIndex(x0 + 1, y0,     z0 + 1, params)];
        const float f011 = maps[base + gridIndex(x0,     y0 + 1, z0 + 1, params)];
        const float f111 = maps[base + gridIndex(x0 + 1, y0 + 1, z0 + 1, params)];

        float value =
            f000 * inverse.x  * inverse.y  * inverse.z +
            f100 * fraction.x * inverse.y  * inverse.z +
            f010 * inverse.x  * fraction.y * inverse.z +
            f110 * fraction.x * fraction.y * inverse.z +
            f001 * inverse.x  * inverse.y  * fraction.z +
            f101 * fraction.x * inverse.y  * fraction.z +
            f011 * inverse.x  * fraction.y * fraction.z +
            f111 * fraction.x * fraction.y * fraction.z;

        float3 gradient = float3(
            f000 * -inverse.y * inverse.z + f100 * inverse.y * inverse.z
                + f010 * -fraction.y * inverse.z + f110 * fraction.y * inverse.z
                + f001 * -inverse.y * fraction.z + f101 * inverse.y * fraction.z
                + f011 * -fraction.y * fraction.z + f111 * fraction.y * fraction.z,
            f000 * inverse.x * -inverse.z + f100 * fraction.x * -inverse.z
                + f010 * inverse.x * inverse.z + f110 * fraction.x * inverse.z
                + f001 * inverse.x * -fraction.z + f101 * fraction.x * -fraction.z
                + f011 * inverse.x * fraction.z + f111 * fraction.x * fraction.z,
            f000 * inverse.x * inverse.y * -1.0f + f100 * fraction.x * inverse.y * -1.0f
                + f010 * inverse.x * fraction.y * -1.0f + f110 * fraction.x * fraction.y * -1.0f
                + f001 * inverse.x * inverse.y + f101 * fraction.x * inverse.y
                + f011 * inverse.x * fraction.y + f111 * fraction.x * fraction.y
        );

        if (value > 0.0f) {
            const float curl = params.curlV / (params.curlV + value);
            value *= curl;
            gradient *= curl * curl;
        }

        energy += value + params.slope * (miss.x + miss.y + miss.z) * params.spacing;
        totalGradient += gradient / params.spacing + params.slope * float3(region);
    }

    scores[pose] = packed_float4(totalGradient, energy);
}

struct SearchParams {
    uint lanes;
    uint steps;
    uint seed;
    uint localSteps;
    uint torsionCount;
    uint pairCount;
    float centerX;
    float centerY;
    float centerZ;
    float padding0;
    float spanX;
    float spanY;
    float spanZ;
    float translationMutation;
    float rotationMutation;
    float temperature;
    float gradientStep;
    float padding1;
};

struct TorsionData {
    uint parent;
    uint child;
    uint maskLow;
    uint maskHigh;
};

struct PairData {
    uint a;
    uint b;
    uint typeA;
    uint typeB;
};

inline uint randomUInt(thread uint& state) {
    state ^= state << 13;
    state ^= state >> 17;
    state ^= state << 5;
    return state;
}

inline float randomUnit(thread uint& state) {
    return float(randomUInt(state) & 0x00ffffffu) / float(0x01000000u);
}

inline float randomNormal(thread uint& state) {
    const float u1 = max(randomUnit(state), 1e-7f);
    const float u2 = randomUnit(state);
    return sqrt(-2.0f * log(u1)) * cos(2.0f * M_PI_F * u2);
}

inline float3 randomInsideSphere(thread uint& state) {
    for (;;) {
        const float3 value = float3(randomUnit(state), randomUnit(state), randomUnit(state)) * 2.0f - 1.0f;
        if (dot(value, value) < 1.0f) return value;
    }
}

inline float4 quaternionMultiply(float4 a, float4 b) {
    return float4(
        a.w * b.xyz + b.w * a.xyz + cross(a.xyz, b.xyz),
        a.w * b.w - dot(a.xyz, b.xyz)
    );
}

inline float3 quaternionRotate(float4 q, float3 v) {
    return v + 2.0f * cross(q.xyz, cross(q.xyz, v) + q.w * v);
}

inline float scoreRigidPose(
    device const float* maps,
    device const packed_float4* atoms,
    float3 position,
    float4 orientation,
    constant GridParams& params,
    thread float3& totalGradient
) {
    const uint mapStride = params.dimX * params.dimY * params.dimZ;
    float energy = 0.0f;
    totalGradient = 0.0f;
    for (uint atomIndex = 0; atomIndex < params.atomCount; ++atomIndex) {
        const packed_float4 atom = atoms[atomIndex];
        const float3 location = quaternionRotate(orientation, atom.xyz) + position;
        const uint mapIndex = min(uint(atom.w), params.mapCount - 1);
        float3 s = (location - float3(params.originX, params.originY, params.originZ)) / params.spacing;
        float3 miss = 0.0f;
        int3 region = 0;
        if (s.x < 0.0f) { miss.x = -s.x; region.x = -1; s.x = 0.0f; }
        else if (s.x >= float(params.dimX - 1)) { miss.x = s.x - float(params.dimX - 1); region.x = 1; s.x = float(params.dimX - 1); }
        if (s.y < 0.0f) { miss.y = -s.y; region.y = -1; s.y = 0.0f; }
        else if (s.y >= float(params.dimY - 1)) { miss.y = s.y - float(params.dimY - 1); region.y = 1; s.y = float(params.dimY - 1); }
        if (s.z < 0.0f) { miss.z = -s.z; region.z = -1; s.z = 0.0f; }
        else if (s.z >= float(params.dimZ - 1)) { miss.z = s.z - float(params.dimZ - 1); region.z = 1; s.z = float(params.dimZ - 1); }
        const uint x0 = min(uint(s.x), params.dimX - 2);
        const uint y0 = min(uint(s.y), params.dimY - 2);
        const uint z0 = min(uint(s.z), params.dimZ - 2);
        const float3 f = clamp(s - float3(x0, y0, z0), 0.0f, 1.0f);
        const float3 m = 1.0f - f;
        const uint base = mapIndex * mapStride;
        const float f000 = maps[base + gridIndex(x0, y0, z0, params)];
        const float f100 = maps[base + gridIndex(x0 + 1, y0, z0, params)];
        const float f010 = maps[base + gridIndex(x0, y0 + 1, z0, params)];
        const float f110 = maps[base + gridIndex(x0 + 1, y0 + 1, z0, params)];
        const float f001 = maps[base + gridIndex(x0, y0, z0 + 1, params)];
        const float f101 = maps[base + gridIndex(x0 + 1, y0, z0 + 1, params)];
        const float f011 = maps[base + gridIndex(x0, y0 + 1, z0 + 1, params)];
        const float f111 = maps[base + gridIndex(x0 + 1, y0 + 1, z0 + 1, params)];
        float value = f000*m.x*m.y*m.z + f100*f.x*m.y*m.z + f010*m.x*f.y*m.z + f110*f.x*f.y*m.z
            + f001*m.x*m.y*f.z + f101*f.x*m.y*f.z + f011*m.x*f.y*f.z + f111*f.x*f.y*f.z;
        float3 gradient = float3(
            f000*-m.y*m.z + f100*m.y*m.z + f010*-f.y*m.z + f110*f.y*m.z + f001*-m.y*f.z + f101*m.y*f.z + f011*-f.y*f.z + f111*f.y*f.z,
            f000*m.x*-m.z + f100*f.x*-m.z + f010*m.x*m.z + f110*f.x*m.z + f001*m.x*-f.z + f101*f.x*-f.z + f011*m.x*f.z + f111*f.x*f.z,
            f000*m.x*m.y*-1.0f + f100*f.x*m.y*-1.0f + f010*m.x*f.y*-1.0f + f110*f.x*f.y*-1.0f + f001*m.x*m.y + f101*f.x*m.y + f011*m.x*f.y + f111*f.x*f.y
        );
        if (value > 0.0f) {
            const float curl = params.curlV / (params.curlV + value);
            value *= curl;
            gradient *= curl * curl;
        }
        energy += value + params.slope * (miss.x + miss.y + miss.z) * params.spacing;
        totalGradient += gradient / params.spacing + params.slope * float3(region);
    }
    return energy;
}

kernel void rigidDock(
    device const float* maps [[buffer(0)]],
    device const packed_float4* atoms [[buffer(1)]],
    device packed_float4* resultPoses [[buffer(2)]],
    device packed_float4* resultOrientations [[buffer(3)]],
    constant GridParams& grid [[buffer(4)]],
    constant SearchParams& search [[buffer(5)]],
    uint lane [[thread_position_in_grid]]
) {
    if (lane >= search.lanes) return;
    uint rng = search.seed ^ (lane * 747796405u + 2891336453u);
    float3 position = float3(search.centerX, search.centerY, search.centerZ)
        + (float3(randomUnit(rng), randomUnit(rng), randomUnit(rng)) - 0.5f)
        * float3(search.spanX, search.spanY, search.spanZ);
    float4 orientation = normalize(float4(
        randomUnit(rng) * 2.0f - 1.0f,
        randomUnit(rng) * 2.0f - 1.0f,
        randomUnit(rng) * 2.0f - 1.0f,
        randomUnit(rng) * 2.0f - 1.0f
    ));
    float3 gradient;
    float energy = scoreRigidPose(maps, atoms, position, orientation, grid, gradient);
    float bestEnergy = energy;
    float3 bestPosition = position;
    float4 bestOrientation = orientation;

    for (uint step = 0; step < search.steps; ++step) {
        float3 candidatePosition = position + (float3(randomUnit(rng), randomUnit(rng), randomUnit(rng)) - 0.5f) * 2.0f * search.translationMutation;
        float3 axis = normalize(float3(randomUnit(rng), randomUnit(rng), randomUnit(rng)) * 2.0f - 1.0f + 1e-6f);
        const float angle = (randomUnit(rng) * 2.0f - 1.0f) * search.rotationMutation;
        const float4 delta = float4(axis * sin(angle * 0.5f), cos(angle * 0.5f));
        float4 candidateOrientation = normalize(quaternionMultiply(delta, orientation));
        float candidateEnergy = scoreRigidPose(maps, atoms, candidatePosition, candidateOrientation, grid, gradient);
        for (uint local = 0; local < search.localSteps; ++local) {
            const float gradientLength = length(gradient);
            if (gradientLength < 1e-5f) break;
            candidatePosition -= gradient / gradientLength * min(search.gradientStep, gradientLength * 0.002f);
            candidateEnergy = scoreRigidPose(maps, atoms, candidatePosition, candidateOrientation, grid, gradient);
        }
        const bool accept = candidateEnergy < energy || randomUnit(rng) < exp((energy - candidateEnergy) / search.temperature);
        if (accept) {
            position = candidatePosition;
            orientation = candidateOrientation;
            energy = candidateEnergy;
            if (energy < bestEnergy) {
                bestEnergy = energy;
                bestPosition = position;
                bestOrientation = orientation;
            }
        }
    }
    resultPoses[lane] = packed_float4(bestPosition, bestEnergy);
    resultOrientations[lane] = packed_float4(bestOrientation);
}

inline bool torsionContains(const TorsionData torsion, uint atom) {
    return atom < 32 ? (torsion.maskLow & (1u << atom)) != 0 : (torsion.maskHigh & (1u << (atom - 32))) != 0;
}

inline float3 rotateAroundAxis(float3 point, float3 axisPoint, float3 axis, float angle) {
    const float3 relative = point - axisPoint;
    const float c = cos(angle);
    const float s = sin(angle);
    return axisPoint + relative * c + cross(axis, relative) * s + axis * dot(axis, relative) * (1.0f - c);
}

inline void buildFlexibleCoordinates(
    device const packed_float4* atoms,
    device const TorsionData* torsions,
    thread const float* angles,
    uint atomCount,
    uint torsionCount,
    thread packed_float4* coordinates
) {
    for (uint atom = 0; atom < atomCount; ++atom) coordinates[atom] = atoms[atom];
    for (uint torsionIndex = 0; torsionIndex < torsionCount; ++torsionIndex) {
        const TorsionData torsion = torsions[torsionIndex];
        const float3 axisPoint = coordinates[torsion.parent].xyz;
        const float3 axis = normalize(coordinates[torsion.child].xyz - axisPoint);
        for (uint atom = 0; atom < atomCount; ++atom) {
            if (torsionContains(torsion, atom)) {
                coordinates[atom].xyz = rotateAroundAxis(coordinates[atom].xyz, axisPoint, axis, angles[torsionIndex]);
            }
        }
    }
}

constant float xsRadiiMetal[16] = {1.9f,1.9f,1.8f,1.8f,1.8f,1.8f,1.7f,1.7f,1.7f,1.7f,2.0f,2.1f,1.5f,1.8f,2.0f,2.2f};
inline float xsRadius(uint type) { return xsRadiiMetal[min(type, 15u)]; }

inline bool xsHydrophobic(uint type) {
    return type == 0 || type == 12 || type == 13 || type == 14 || type == 15;
}

inline bool xsAcceptor(uint type) { return type == 4 || type == 5 || type == 8 || type == 9; }
inline bool xsDonor(uint type) { return type == 3 || type == 5 || type == 7 || type == 9 || type == 18; }
inline bool xsHydrogenBond(uint a, uint b) { return (xsDonor(a) && xsAcceptor(b)) || (xsDonor(b) && xsAcceptor(a)); }

inline float slopeStep(float bad, float good, float value) {
    if (bad < good) {
        if (value <= bad) return 0.0f;
        if (value >= good) return 1.0f;
    } else {
        if (value >= bad) return 0.0f;
        if (value <= good) return 1.0f;
    }
    return (value - bad) / (good - bad);
}

inline float intramolecularEnergy(thread const packed_float4* coordinates, device const PairData* pairs, uint pairCount) {
    float energy = 0.0f;
    for (uint index = 0; index < pairCount; ++index) {
        const PairData pair = pairs[index];
        const float distance = length(coordinates[pair.a].xyz - coordinates[pair.b].xyz);
        if (distance >= 8.0f) continue;
        const float d = distance - xsRadius(pair.typeA) - xsRadius(pair.typeB);
        const float gauss1 = exp(-pow(d / 0.5f, 2.0f));
        const float gauss2 = exp(-pow((d - 3.0f) / 2.0f, 2.0f));
        const float repulsion = d < 0.0f ? d * d : 0.0f;
        const float hydrophobic = xsHydrophobic(pair.typeA) && xsHydrophobic(pair.typeB) ? slopeStep(1.5f, 0.5f, d) : 0.0f;
        const float hydrogen = xsHydrogenBond(pair.typeA, pair.typeB) ? slopeStep(0.0f, -0.7f, d) : 0.0f;
        energy += -0.035579f * gauss1 - 0.005156f * gauss2 + 0.840245f * repulsion
            - 0.035069f * hydrophobic - 0.587439f * hydrogen;
    }
    return energy;
}

inline float slopeDerivative(float bad, float good, float value) {
    const float low = min(bad, good), high = max(bad, good);
    return value > low && value < high ? 1.0f / (good - bad) : 0.0f;
}

inline float intramolecularEnergyGradient(
    thread const packed_float4* coordinates,
    device const PairData* pairs,
    uint pairCount,
    thread float3* gradients
) {
    float energy = 0.0f;
    for (uint index = 0; index < pairCount; ++index) {
        const PairData pair = pairs[index];
        const float3 delta = coordinates[pair.a].xyz - coordinates[pair.b].xyz;
        const float distance = length(delta);
        if (distance >= 8.0f || distance < 1e-6f) continue;
        const float d = distance - xsRadius(pair.typeA) - xsRadius(pair.typeB);
        const float gauss1 = exp(-pow(d / 0.5f, 2.0f));
        const float gauss2 = exp(-pow((d - 3.0f) / 2.0f, 2.0f));
        const bool hydroPair = xsHydrophobic(pair.typeA) && xsHydrophobic(pair.typeB);
        const bool hydrogenPair = xsHydrogenBond(pair.typeA, pair.typeB);
        const float repulsion = d < 0.0f ? d * d : 0.0f;
        const float hydrophobic = hydroPair ? slopeStep(1.5f, 0.5f, d) : 0.0f;
        const float hydrogen = hydrogenPair ? slopeStep(0.0f, -0.7f, d) : 0.0f;
        energy += -0.035579f * gauss1 - 0.005156f * gauss2 + 0.840245f * repulsion
            - 0.035069f * hydrophobic - 0.587439f * hydrogen;
        float derivative = -0.035579f * gauss1 * (-8.0f * d)
            - 0.005156f * gauss2 * (-(d - 3.0f) * 0.5f);
        if (d < 0.0f) derivative += 0.840245f * 2.0f * d;
        if (hydroPair) derivative += -0.035069f * slopeDerivative(1.5f, 0.5f, d);
        if (hydrogenPair) derivative += -0.587439f * slopeDerivative(0.0f, -0.7f, d);
        const float3 cartesian = derivative * delta / distance;
        gradients[pair.a] += cartesian;
        gradients[pair.b] -= cartesian;
    }
    return energy;
}

inline float scoreFlexiblePose(
    device const float* maps,
    thread const packed_float4* coordinates,
    device const PairData* pairs,
    float3 position,
    float4 orientation,
    constant GridParams& params,
    constant SearchParams& search,
    thread float3& totalGradient
) {
    const uint mapStride = params.dimX * params.dimY * params.dimZ;
    float energy = intramolecularEnergy(coordinates, pairs, search.pairCount);
    totalGradient = 0.0f;
    for (uint atomIndex = 0; atomIndex < params.atomCount; ++atomIndex) {
        const packed_float4 atom = coordinates[atomIndex];
        const float3 location = quaternionRotate(orientation, atom.xyz) + position;
        const uint mapIndex = min(uint(atom.w), params.mapCount - 1);
        float3 s = (location - float3(params.originX, params.originY, params.originZ)) / params.spacing;
        float3 miss = 0.0f;
        int3 region = 0;
        if (s.x < 0.0f) { miss.x = -s.x; region.x = -1; s.x = 0.0f; } else if (s.x >= float(params.dimX - 1)) { miss.x = s.x - float(params.dimX - 1); region.x = 1; s.x = float(params.dimX - 1); }
        if (s.y < 0.0f) { miss.y = -s.y; region.y = -1; s.y = 0.0f; } else if (s.y >= float(params.dimY - 1)) { miss.y = s.y - float(params.dimY - 1); region.y = 1; s.y = float(params.dimY - 1); }
        if (s.z < 0.0f) { miss.z = -s.z; region.z = -1; s.z = 0.0f; } else if (s.z >= float(params.dimZ - 1)) { miss.z = s.z - float(params.dimZ - 1); region.z = 1; s.z = float(params.dimZ - 1); }
        const uint x0 = min(uint(s.x), params.dimX - 2), y0 = min(uint(s.y), params.dimY - 2), z0 = min(uint(s.z), params.dimZ - 2);
        const float3 f = clamp(s - float3(x0, y0, z0), 0.0f, 1.0f), m = 1.0f - f;
        const uint base = mapIndex * mapStride;
        const float f000=maps[base+gridIndex(x0,y0,z0,params)], f100=maps[base+gridIndex(x0+1,y0,z0,params)];
        const float f010=maps[base+gridIndex(x0,y0+1,z0,params)], f110=maps[base+gridIndex(x0+1,y0+1,z0,params)];
        const float f001=maps[base+gridIndex(x0,y0,z0+1,params)], f101=maps[base+gridIndex(x0+1,y0,z0+1,params)];
        const float f011=maps[base+gridIndex(x0,y0+1,z0+1,params)], f111=maps[base+gridIndex(x0+1,y0+1,z0+1,params)];
        float value=f000*m.x*m.y*m.z+f100*f.x*m.y*m.z+f010*m.x*f.y*m.z+f110*f.x*f.y*m.z+f001*m.x*m.y*f.z+f101*f.x*m.y*f.z+f011*m.x*f.y*f.z+f111*f.x*f.y*f.z;
        float3 gradient=float3(
            f000*-m.y*m.z+f100*m.y*m.z+f010*-f.y*m.z+f110*f.y*m.z+f001*-m.y*f.z+f101*m.y*f.z+f011*-f.y*f.z+f111*f.y*f.z,
            f000*m.x*-m.z+f100*f.x*-m.z+f010*m.x*m.z+f110*f.x*m.z+f001*m.x*-f.z+f101*f.x*-f.z+f011*m.x*f.z+f111*f.x*f.z,
            f000*m.x*m.y*-1.0f+f100*f.x*m.y*-1.0f+f010*m.x*f.y*-1.0f+f110*f.x*f.y*-1.0f+f001*m.x*m.y+f101*f.x*m.y+f011*m.x*f.y+f111*f.x*f.y);
        if(value>0.0f){const float curl=params.curlV/(params.curlV+value);value*=curl;gradient*=curl*curl;}
        energy+=value+params.slope*(miss.x+miss.y+miss.z)*params.spacing;
        totalGradient+=gradient/params.spacing+params.slope*float3(region);
    }
    return energy;
}

inline void flexibleDofGradient(
    device const float* maps,
    thread const packed_float4* coordinates,
    device const TorsionData* torsions,
    device const PairData* pairs,
    float3 position,
    float4 orientation,
    constant GridParams& params,
    constant SearchParams& search,
    thread float* dofGradient
) {
    float3 atomGradient[64];
    for (uint atom = 0; atom < params.atomCount; ++atom) {
        atomGradient[atom] = 0.0f;
    }
    intramolecularEnergyGradient(coordinates, pairs, search.pairCount, atomGradient);
    const uint mapStride = params.dimX * params.dimY * params.dimZ;
    for (uint atomIndex = 0; atomIndex < params.atomCount; ++atomIndex) {
        const packed_float4 atom = coordinates[atomIndex];
        const float3 location = quaternionRotate(orientation, atom.xyz) + position;
        const uint mapIndex = min(uint(atom.w), params.mapCount - 1);
        float3 s = (location - float3(params.originX, params.originY, params.originZ)) / params.spacing;
        int3 region = 0;
        if (s.x < 0.0f) { region.x = -1; s.x = 0.0f; } else if (s.x >= float(params.dimX - 1)) { region.x = 1; s.x = float(params.dimX - 1); }
        if (s.y < 0.0f) { region.y = -1; s.y = 0.0f; } else if (s.y >= float(params.dimY - 1)) { region.y = 1; s.y = float(params.dimY - 1); }
        if (s.z < 0.0f) { region.z = -1; s.z = 0.0f; } else if (s.z >= float(params.dimZ - 1)) { region.z = 1; s.z = float(params.dimZ - 1); }
        const uint x0 = min(uint(s.x), params.dimX - 2), y0 = min(uint(s.y), params.dimY - 2), z0 = min(uint(s.z), params.dimZ - 2);
        const float3 f = clamp(s - float3(x0, y0, z0), 0.0f, 1.0f), m = 1.0f - f;
        const uint base = mapIndex * mapStride;
        const float f000=maps[base+gridIndex(x0,y0,z0,params)], f100=maps[base+gridIndex(x0+1,y0,z0,params)];
        const float f010=maps[base+gridIndex(x0,y0+1,z0,params)], f110=maps[base+gridIndex(x0+1,y0+1,z0,params)];
        const float f001=maps[base+gridIndex(x0,y0,z0+1,params)], f101=maps[base+gridIndex(x0+1,y0,z0+1,params)];
        const float f011=maps[base+gridIndex(x0,y0+1,z0+1,params)], f111=maps[base+gridIndex(x0+1,y0+1,z0+1,params)];
        float value=f000*m.x*m.y*m.z+f100*f.x*m.y*m.z+f010*m.x*f.y*m.z+f110*f.x*f.y*m.z+f001*m.x*m.y*f.z+f101*f.x*m.y*f.z+f011*m.x*f.y*f.z+f111*f.x*f.y*f.z;
        float3 gridGradient=float3(
            f000*-m.y*m.z+f100*m.y*m.z+f010*-f.y*m.z+f110*f.y*m.z+f001*-m.y*f.z+f101*m.y*f.z+f011*-f.y*f.z+f111*f.y*f.z,
            f000*m.x*-m.z+f100*f.x*-m.z+f010*m.x*m.z+f110*f.x*m.z+f001*m.x*-f.z+f101*f.x*-f.z+f011*m.x*f.z+f111*f.x*f.z,
            f000*m.x*m.y*-1.0f+f100*f.x*m.y*-1.0f+f010*m.x*f.y*-1.0f+f110*f.x*f.y*-1.0f+f001*m.x*m.y+f101*f.x*m.y+f011*m.x*f.y+f111*f.x*f.y);
        if(value>0.0f){const float curl=params.curlV/(params.curlV+value);gridGradient*=curl*curl;}
        atomGradient[atomIndex] = quaternionRotate(orientation, atomGradient[atomIndex])
            + gridGradient / params.spacing + params.slope * float3(region);
    }
    const uint dimensions = 6 + search.torsionCount;
    for (uint i = 0; i < dimensions; ++i) dofGradient[i] = 0.0f;
    for (uint atom = 0; atom < params.atomCount; ++atom) {
        dofGradient[0] += atomGradient[atom].x;
        dofGradient[1] += atomGradient[atom].y;
        dofGradient[2] += atomGradient[atom].z;
        const float3 relative = quaternionRotate(orientation, coordinates[atom].xyz);
        const float3 torque = cross(relative, atomGradient[atom]);
        dofGradient[3] += torque.x;
        dofGradient[4] += torque.y;
        dofGradient[5] += torque.z;
    }
    for (uint torsionIndex = 0; torsionIndex < search.torsionCount; ++torsionIndex) {
        const TorsionData torsion = torsions[torsionIndex];
        const float3 axisPoint = coordinates[torsion.parent].xyz;
        const float3 axis = normalize(coordinates[torsion.child].xyz - axisPoint);
        float derivative = 0.0f;
        for (uint atom = 0; atom < params.atomCount; ++atom) {
            if (!torsionContains(torsion, atom)) continue;
            const float3 localVelocity = cross(axis, coordinates[atom].xyz - axisPoint);
            derivative += dot(atomGradient[atom], quaternionRotate(orientation, localVelocity));
        }
        dofGradient[6 + torsionIndex] = derivative;
    }
}

inline void applyDofStep(
    float3 position,
    float4 orientation,
    thread const float* angles,
    thread const float* step,
    uint torsionCount,
    thread float3& resultPosition,
    thread float4& resultOrientation,
    thread float* resultAngles
) {
    resultPosition = position + float3(step[0], step[1], step[2]);
    const float3 rotation = float3(step[3], step[4], step[5]);
    const float rotationAngle = length(rotation);
    resultOrientation = orientation;
    if (rotationAngle > 1e-8f) {
        const float halfAngle = 0.5f * rotationAngle;
        resultOrientation = normalize(quaternionMultiply(float4(rotation / rotationAngle * sin(halfAngle), cos(halfAngle)), orientation));
    }
    for (uint i = 0; i < 8; ++i) resultAngles[i] = angles[i] + (i < torsionCount ? step[6 + i] : 0.0f);
}

inline uint symmetricIndex(uint row, uint column) {
    const uint high = max(row, column), low = min(row, column);
    return high * (high + 1u) / 2u + low;
}

inline float bfgsRefine(
    device const float* maps,
    device const packed_float4* atoms,
    device const TorsionData* torsions,
    device const PairData* pairs,
    constant GridParams& grid,
    constant SearchParams& search,
    thread float3& position,
    thread float4& orientation,
    thread float* angles
) {
    const uint n = 6 + search.torsionCount;
    float inverseHessian[105];
    float gradient[14], nextGradient[14], direction[14], stepVector[14], y[14], hy[14];
    for (uint row = 0; row < n; ++row) for (uint column = 0; column <= row; ++column)
        inverseHessian[symmetricIndex(row, column)] = row == column ? 1.0f : 0.0f;
    packed_float4 coordinates[64];
    const float3 originalPosition = position;
    const float4 originalOrientation = orientation;
    float originalAngles[8]; for (uint i = 0; i < 8; ++i) originalAngles[i] = angles[i];
    buildFlexibleCoordinates(atoms,torsions,angles,grid.atomCount,search.torsionCount,coordinates);
    float3 ignored;
    float energy = scoreFlexiblePose(maps,coordinates,pairs,position,orientation,grid,search,ignored);
    const float originalEnergy = energy;
    flexibleDofGradient(maps,coordinates,torsions,pairs,position,orientation,grid,search,gradient);
    for (uint iteration = 0; iteration < search.localSteps; ++iteration) {
        float gradientNorm = 0.0f;
        for (uint i = 0; i < n; ++i) gradientNorm += gradient[i] * gradient[i];
        if (gradientNorm < 1e-8f) break;
        for (uint row = 0; row < n; ++row) {
            direction[row] = 0.0f;
            for (uint column = 0; column < n; ++column) direction[row] -= inverseHessian[symmetricIndex(row, column)] * gradient[column];
        }
        float directionGradient = 0.0f;
        for (uint i = 0; i < n; ++i) directionGradient += direction[i] * gradient[i];
        float3 trialPosition = position;
        float4 trialOrientation = orientation;
        float trialAngles[8];
        float trialEnergy = energy;
        float alpha = 1.0f;
        for (uint lineSearch = 0; lineSearch < 10; ++lineSearch) {
            for (uint i = 0; i < n; ++i) stepVector[i] = direction[i] * alpha;
            applyDofStep(position,orientation,angles,stepVector,search.torsionCount,trialPosition,trialOrientation,trialAngles);
            buildFlexibleCoordinates(atoms,torsions,trialAngles,grid.atomCount,search.torsionCount,coordinates);
            trialEnergy = scoreFlexiblePose(maps,coordinates,pairs,trialPosition,trialOrientation,grid,search,ignored);
            if (trialEnergy - energy < 0.0001f * alpha * directionGradient) break;
            alpha *= 0.5f;
        }
        flexibleDofGradient(maps,coordinates,torsions,pairs,trialPosition,trialOrientation,grid,search,nextGradient);
        float ys = 0.0f;
        for (uint i = 0; i < n; ++i) { y[i] = nextGradient[i] - gradient[i]; ys += y[i] * stepVector[i]; }
        if (iteration == 0) {
            float yy = 0.0f;
            for (uint i = 0; i < n; ++i) yy += y[i] * y[i];
            if (abs(yy) > 1e-7f) {
                const float diagonal = ys / yy;
                for (uint row = 0; row < n; ++row) for (uint column = 0; column <= row; ++column)
                    inverseHessian[symmetricIndex(row, column)] = row == column ? diagonal : 0.0f;
            }
        }
        if (ys > 1e-7f) {
            float yHy = 0.0f;
            for (uint row = 0; row < n; ++row) {
                hy[row] = 0.0f;
                for (uint column = 0; column < n; ++column) hy[row] += inverseHessian[symmetricIndex(row, column)] * y[column];
                yHy += y[row] * hy[row];
            }
            const float factor = (1.0f + yHy / ys) / ys;
            for (uint row = 0; row < n; ++row) for (uint column = 0; column <= row; ++column)
                inverseHessian[symmetricIndex(row, column)] += factor * stepVector[row] * stepVector[column]
                    - (stepVector[row] * hy[column] + hy[row] * stepVector[column]) / ys;
        } else {
            for (uint row = 0; row < n; ++row) for (uint column = 0; column <= row; ++column)
                inverseHessian[symmetricIndex(row, column)] = row == column ? 1.0f : 0.0f;
        }
        position = trialPosition; orientation = trialOrientation; energy = trialEnergy;
        for (uint i = 0; i < 8; ++i) angles[i] = trialAngles[i];
        for (uint i = 0; i < n; ++i) gradient[i] = nextGradient[i];
    }
    if (energy > originalEnergy) {
        position = originalPosition; orientation = originalOrientation;
        for (uint i = 0; i < 8; ++i) angles[i] = originalAngles[i];
        buildFlexibleCoordinates(atoms,torsions,angles,grid.atomCount,search.torsionCount,coordinates);
        energy = originalEnergy;
    }
    return energy;
}

kernel void flexibleDock(
    device const float* maps [[buffer(0)]],
    device const packed_float4* atoms [[buffer(1)]],
    device packed_float4* resultPoses [[buffer(2)]],
    device packed_float4* resultOrientations [[buffer(3)]],
    constant GridParams& grid [[buffer(4)]],
    constant SearchParams& search [[buffer(5)]],
    device const TorsionData* torsions [[buffer(6)]],
    device const PairData* pairs [[buffer(7)]],
    device packed_float4* resultAngles [[buffer(8)]],
    uint lane [[thread_position_in_grid]]
) {
    if (lane >= search.lanes || grid.atomCount > 64 || search.torsionCount > 8) return;
    uint rng = search.seed ^ (lane * 747796405u + 2891336453u);
    float3 position=float3(search.centerX,search.centerY,search.centerZ)+(float3(randomUnit(rng),randomUnit(rng),randomUnit(rng))-0.5f)*float3(search.spanX,search.spanY,search.spanZ);
    float4 orientation=normalize(float4(randomNormal(rng),randomNormal(rng),randomNormal(rng),randomNormal(rng)));
    float angles[8]; for(uint i=0;i<8;++i)angles[i]=i<search.torsionCount ? (randomUnit(rng)*2.0f-1.0f)*M_PI_F : 0.0f;
    packed_float4 coordinates[64];
    buildFlexibleCoordinates(atoms,torsions,angles,grid.atomCount,search.torsionCount,coordinates);
    float3 gradient;
    float energy=scoreFlexiblePose(maps,coordinates,pairs,position,orientation,grid,search,gradient);
    float bestEnergy=energy; float3 bestPosition=position; float4 bestOrientation=orientation; float bestAngles[8]={0,0,0,0,0,0,0,0};
    for(uint step=0;step<search.steps;++step){
        float3 candidatePosition=position;
        float4 candidateOrientation=orientation;
        float candidateAngles[8]; for(uint i=0;i<8;++i)candidateAngles[i]=angles[i];
        const uint changed=randomUInt(rng)%(2u+search.torsionCount);
        if(changed==0u) candidatePosition += search.translationMutation * randomInsideSphere(rng);
        else if(changed==1u) {
            const float3 rotation = search.rotationMutation * randomInsideSphere(rng);
            const float rotationAngle=length(rotation);
            if(rotationAngle>1e-8f) candidateOrientation=normalize(quaternionMultiply(float4(rotation/rotationAngle*sin(rotationAngle*0.5f),cos(rotationAngle*0.5f)),orientation));
        } else candidateAngles[changed-2u]=(randomUnit(rng)*2.0f-1.0f)*M_PI_F;
        float candidateEnergy=bfgsRefine(maps,atoms,torsions,pairs,grid,search,candidatePosition,candidateOrientation,candidateAngles);
        const bool accept=step==0u||candidateEnergy<energy||randomUnit(rng)<exp((energy-candidateEnergy)/search.temperature);
        if(accept){position=candidatePosition;orientation=candidateOrientation;energy=candidateEnergy;for(uint i=0;i<8;++i)angles[i]=candidateAngles[i];if(energy<bestEnergy){bestEnergy=energy;bestPosition=position;bestOrientation=orientation;for(uint i=0;i<8;++i)bestAngles[i]=angles[i];}}
    }
    resultPoses[lane]=packed_float4(bestPosition,bestEnergy);resultOrientations[lane]=packed_float4(bestOrientation);
    resultAngles[lane*2]=packed_float4(bestAngles[0],bestAngles[1],bestAngles[2],bestAngles[3]);
    resultAngles[lane*2+1]=packed_float4(bestAngles[4],bestAngles[5],bestAngles[6],bestAngles[7]);
}
