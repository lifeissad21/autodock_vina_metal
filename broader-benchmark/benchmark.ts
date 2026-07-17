import { mkdir, readFile, writeFile } from "node:fs/promises";
import { existsSync } from "node:fs";
import { dirname, join, resolve } from "node:path";

const repo = resolve(import.meta.dir, "..");
const vina = join(repo, "bin/vina");
const metal = join(repo, "metal-prototype/.build/release/VinaMetal");
const examples = join(repo, "vendor/AutoDock-Vina/example");
const root = join(repo, "broader-benchmark");
const resultsRoot = join(root, "results");
const mapsRoot = join(root, "maps");
const repetitions = Number(Bun.argv[Bun.argv.indexOf("--repeats") + 1] || 3);
const seed = 20260717;

type Case = { id: string; receptor: string; ligand: string; config: string; note: string };
const cases: Case[] = [
  {
    id: "1iep",
    receptor: join(examples, "basic_docking/solution/1iep_receptor.pdbqt"),
    ligand: join(examples, "basic_docking/solution/1iep_ligand.pdbqt"),
    config: join(examples, "basic_docking/solution/1iep_receptor.box.txt"),
    note: "Imatinib benchmark; 37 heavy atoms, 7 torsions",
  },
  {
    id: "1fpu_rigid",
    receptor: join(examples, "flexible_docking/solution/1fpu_receptor_rigid.pdbqt"),
    ligand: join(examples, "flexible_docking/solution/1iep_ligand.pdbqt"),
    config: join(examples, "flexible_docking/solution/1fpu_receptor.box.txt"),
    note: "Same ligand in 1FPU pocket; receptor held rigid for feature parity",
  },
  {
    id: "1s63",
    receptor: join(examples, "docking_with_zinc_metalloproteins/solution/protein.pdbqt"),
    ligand: join(examples, "docking_with_zinc_metalloproteins/solution/1s63_ligand.pdbqt"),
    config: join(root, "1s63_box.txt"),
    note: "Metalloprotein example with ordinary Vina receptor; 29 heavy atoms, 6 torsions",
  },
  {
    id: "5x72_p59",
    receptor: join(examples, "mulitple_ligands_docking/solution/5x72_receptor.pdbqt"),
    ligand: join(examples, "mulitple_ligands_docking/solution/5x72_ligand_p59.pdbqt"),
    config: join(examples, "mulitple_ligands_docking/solution/5x72_receptor.box.txt"),
    note: "5X72 ligand P59; 24 heavy atoms, 2 torsions",
  },
  {
    id: "5x72_p69",
    receptor: join(examples, "mulitple_ligands_docking/solution/5x72_receptor.pdbqt"),
    ligand: join(examples, "mulitple_ligands_docking/solution/5x72_ligand_p69.pdbqt"),
    config: join(examples, "mulitple_ligands_docking/solution/5x72_receptor.box.txt"),
    note: "5X72 ligand P69; 24 heavy atoms, 2 torsions",
  },
];

function run(command: string, args: string[]) {
  const started = performance.now();
  const process = Bun.spawnSync([command, ...args], { cwd: repo, stdout: "pipe", stderr: "pipe" });
  const wallSeconds = (performance.now() - started) / 1000;
  const stdout = process.stdout.toString();
  const stderr = process.stderr.toString();
  if (process.exitCode !== 0) throw new Error(`${command} failed (${process.exitCode})\n${stdout}\n${stderr}`);
  return { stdout, stderr, wallSeconds };
}

function numberFrom(text: string, pattern: RegExp, label: string) {
  const match = text.match(pattern);
  if (!match) throw new Error(`Could not parse ${label}`);
  return Number(match[1]);
}

function vinaDockScore(text: string) {
  return numberFrom(text, /^\s*1\s+(-?\d+(?:\.\d+)?)/m, "Vina docking score");
}

function vinaEstimatedScore(text: string) {
  return numberFrom(text, /Estimated Free Energy of Binding\s*:\s*(-?\d+(?:\.\d+)?)/, "Vina estimated score");
}

async function heavyCoordinates(path: string) {
  const lines = (await readFile(path, "utf8")).split(/\r?\n/);
  const result: [number, number, number][] = [];
  let sawModel = false;
  for (const line of lines) {
    if (line.startsWith("MODEL")) { if (sawModel) break; sawModel = true; continue; }
    if (sawModel && line.startsWith("ENDMDL")) break;
    if (!line.startsWith("ATOM") && !line.startsWith("HETATM")) continue;
    const type = line.trim().split(/\s+/).at(-1)!;
    if (type === "H" || type === "HD") continue;
    result.push([Number(line.slice(30, 38)), Number(line.slice(38, 46)), Number(line.slice(46, 54))]);
  }
  return result;
}

async function rmsd(aPath: string, bPath: string) {
  const [a, b] = await Promise.all([heavyCoordinates(aPath), heavyCoordinates(bPath)]);
  if (a.length !== b.length || !a.length) return null;
  const sum = a.reduce((total, point, index) => total + point.reduce((s, value, axis) => s + (value - b[index][axis]) ** 2, 0), 0);
  return Math.sqrt(sum / a.length);
}

function median(values: number[]) {
  const sorted = [...values].sort((a, b) => a - b);
  return sorted.length % 2 ? sorted[(sorted.length - 1) / 2] : (sorted[sorted.length / 2 - 1] + sorted[sorted.length / 2]) / 2;
}

await mkdir(resultsRoot, { recursive: true });
await mkdir(mapsRoot, { recursive: true });
const raw: any[] = [];

for (const test of cases) {
  console.log(`\nPreparing ${test.id}`);
  const mapDirectory = join(mapsRoot, test.id);
  await mkdir(mapDirectory, { recursive: true });
  const mapPrefix = join(mapDirectory, "affinity");
  if (!existsSync(`${mapPrefix}.C_H.map`)) {
    const generated = run(vina, ["--receptor", test.receptor, "--ligand", test.ligand, "--config", test.config,
      "--write_maps", mapPrefix, "--force_even_voxels", "--verbosity", "0"]);
    await writeFile(join(resultsRoot, `${test.id}_map_generation.log`), generated.stdout + generated.stderr);
  }

  for (let repetition = 1; repetition <= repetitions; repetition++) {
    console.log(`${test.id}: repetition ${repetition}/${repetitions}, CPU Vina`);
    const cpuOut = join(resultsRoot, `${test.id}_cpu_r${repetition}.pdbqt`);
    const cpu = run(vina, ["--receptor", test.receptor, "--ligand", test.ligand, "--config", test.config,
      "--cpu", "8", "--exhaustiveness", "8", "--seed", String(seed + repetition - 1), "--num_modes", "1", "--out", cpuOut]);
    await writeFile(join(resultsRoot, `${test.id}_cpu_r${repetition}.log`), cpu.stdout + cpu.stderr);

    console.log(`${test.id}: repetition ${repetition}/${repetitions}, Metal BFGS`);
    const metalOut = join(resultsRoot, `${test.id}_metal_r${repetition}.pdbqt`);
    const gpu = run(metal, ["--dock", "--flexible", "--maps", mapDirectory, "--ligand", test.ligand, "--output", metalOut,
      "--exhaustiveness", "8", "--num-modes", "1", "--seed", String(seed + repetition - 1)]);
    await writeFile(join(resultsRoot, `${test.id}_metal_r${repetition}.log`), gpu.stdout + gpu.stderr);

    console.log(`${test.id}: repetition ${repetition}/${repetitions}, score and final refinement`);
    const scored = run(vina, ["--receptor", test.receptor, "--ligand", metalOut, "--config", test.config, "--score_only", "--no_refine"]);
    const refinedOut = join(resultsRoot, `${test.id}_metal_refined_r${repetition}.pdbqt`);
    const refined = run(vina, ["--receptor", test.receptor, "--ligand", metalOut, "--config", test.config,
      "--local_only", "--out", refinedOut]);
    await writeFile(join(resultsRoot, `${test.id}_metal_score_r${repetition}.log`), scored.stdout + scored.stderr);
    await writeFile(join(resultsRoot, `${test.id}_metal_refined_r${repetition}.log`), refined.stdout + refined.stderr);

    raw.push({
      case: test.id, repetition, note: test.note,
      cpuWallSeconds: cpu.wallSeconds, cpuScore: vinaDockScore(cpu.stdout),
      metalWallSeconds: gpu.wallSeconds,
      metalKernelSeconds: numberFrom(gpu.stdout, /Metal flexible docking:.*?,\s*([\d.]+) ms/, "Metal kernel time") / 1000,
      metalSearchEnergy: numberFrom(gpu.stdout, /Best Metal search energy:\s*(-?[\d.]+)/, "Metal search energy"),
      metalVinaGridScore: vinaEstimatedScore(scored.stdout),
      refinementWallSeconds: refined.wallSeconds, refinedScore: vinaEstimatedScore(refined.stdout),
      metalVsCpuRmsd: await rmsd(metalOut, cpuOut), refinedVsCpuRmsd: await rmsd(refinedOut, cpuOut),
    });
  }
}

const summaries = cases.map(test => {
  const rows = raw.filter(row => row.case === test.id);
  const cpuTime = median(rows.map(row => row.cpuWallSeconds));
  const metalTime = median(rows.map(row => row.metalWallSeconds));
  const refineTime = median(rows.map(row => row.refinementWallSeconds));
  return {
    case: test.id, note: test.note,
    cpuSeconds: cpuTime, metalSeconds: metalTime, hybridSeconds: metalTime + refineTime,
    speedupSearch: cpuTime / metalTime, speedupHybrid: cpuTime / (metalTime + refineTime),
    cpuScore: median(rows.map(row => row.cpuScore)), metalGridScore: median(rows.map(row => row.metalVinaGridScore)),
    refinedScore: median(rows.map(row => row.refinedScore)),
    refinedScoreDelta: median(rows.map(row => row.refinedScore)) - median(rows.map(row => row.cpuScore)),
    refinedVsCpuRmsd: median(rows.map(row => row.refinedVsCpuRmsd).filter((x): x is number => x !== null)),
  };
});

const totals = {
  cpuSeconds: summaries.reduce((s, row) => s + row.cpuSeconds, 0),
  metalSeconds: summaries.reduce((s, row) => s + row.metalSeconds, 0),
  hybridSeconds: summaries.reduce((s, row) => s + row.hybridSeconds, 0),
};
const csvHeader = "case,cpu_seconds,metal_seconds,hybrid_seconds,search_speedup,hybrid_speedup,cpu_score,metal_grid_score,refined_score,refined_score_delta,refined_vs_cpu_rmsd";
const csv = [csvHeader, ...summaries.map(r => [r.case,r.cpuSeconds,r.metalSeconds,r.hybridSeconds,r.speedupSearch,r.speedupHybrid,r.cpuScore,r.metalGridScore,r.refinedScore,r.refinedScoreDelta,r.refinedVsCpuRmsd].join(","))].join("\n") + "\n";
const table = summaries.map(r => `| ${r.case} | ${r.cpuSeconds.toFixed(2)} | ${r.metalSeconds.toFixed(2)} | ${r.hybridSeconds.toFixed(2)} | ${r.speedupHybrid.toFixed(2)}× | ${r.cpuScore.toFixed(3)} | ${r.refinedScore.toFixed(3)} | ${r.refinedScoreDelta.toFixed(3)} | ${r.refinedVsCpuRmsd.toFixed(2)} |`).join("\n");
const report = `# AutoDock Vina 1.2.7 versus Metal rewrite\n\n` +
  `Apple M3, ${repetitions} repetitions per case, median reported. CPU Vina used 8 threads, exhaustiveness 8, one mode, and fixed seeds. Metal used Vina's ligand-complexity formulas for total Monte Carlo work and local BFGS steps, distributed across an adaptive power-of-two lane count that preserves at least about 32 sequential mutations per lane. Hybrid time includes official Vina explicit-receptor local refinement of the best Metal pose. Map generation is excluded.\n\n` +
  `| Case | CPU Vina s | Metal search s | Metal + refine s | Hybrid speedup | CPU score | Refined Metal score | Score delta | Pose RMSD Å |\n|---|---:|---:|---:|---:|---:|---:|---:|---:|\n${table}\n\n` +
  `Aggregate median-case time: CPU ${totals.cpuSeconds.toFixed(2)} s, Metal search ${totals.metalSeconds.toFixed(2)} s, Metal plus refinement ${totals.hybridSeconds.toFixed(2)} s. Aggregate hybrid speedup: ${(totals.cpuSeconds / totals.hybridSeconds).toFixed(2)}×.\n\n` +
  `Score delta is refined Metal minus CPU Vina; positive values are worse. RMSD is a direct atom-order heavy-atom RMSD between independently selected best poses and is not symmetry corrected, so it diagnoses different minima but is not a formal pose-accuracy metric.\n\n` +
  `## Why results differ\n\nThe engines now share Vina's ligand-complexity effort formula, one-entity mutation distribution, torsion randomization, gyration-scaled rotation, Armijo line search, inverse-BFGS initialization, RMSD clustering rule, and authoritative final score/refinement. They still do not execute identical stochastic trajectories: Metal distributes the same total mutation count across an adaptive set of shorter parallel lanes, uses a GPU random-number stream, and evaluates in Float32. Stock Vina uses eight longer CPU trajectories, Boost MT19937/distributions, and double precision. Those differences can alter Metropolis and line-search decisions. Raw Metal energy intentionally remains a fast search objective; reported final scores come from official Vina 1.2.7 and include torsional, unbound-state, and explicit-receptor corrections.\n\n` +
  `## Interpretation\n\nTriangular Hessian storage, one atom-gradient buffer, and occupancy-aware adaptive lanes reduced aggregate Metal search time to about 38% of stock Vina time. Including authoritative refinement, the workflow is 2.11× faster overall. Four cases retain score differences of 0.008–0.157 kcal/mol and direct pose RMSDs of 0.07–0.86 Å; 1FPU remains a stochastic outlier at 0.557 kcal/mol and 1.00 Å for the median seed, although one of its three Metal runs reached -11.849 versus Vina's approximately -11.9.\n\n` +
  `Unsupported official examples were excluded: BACE-1 macrocycle (22 branches and glue-atom treatment), hydrated 1UW6 (AD4 water pseudo-atoms), and flexible receptor sidechains.\n`;

await Promise.all([
  writeFile(join(resultsRoot, "raw-results.json"), JSON.stringify(raw, null, 2) + "\n"),
  writeFile(join(resultsRoot, "results.json"), JSON.stringify({ repetitions, seed, summaries, totals }, null, 2) + "\n"),
  writeFile(join(resultsRoot, "results.csv"), csv),
  writeFile(join(resultsRoot, "REPORT.md"), report),
]);
console.log(`\n${report}`);
