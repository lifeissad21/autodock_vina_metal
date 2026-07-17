#!/usr/bin/env bun

import { mkdir } from "node:fs/promises";
import { join } from "node:path";

type Run = {
  cpu: number;
  repeat: number;
  seconds: number;
  affinityKcalMol: number | null;
  exitCode: number;
};

const root = import.meta.dir;
const args = new Map<string, string>();
for (let i = 2; i < Bun.argv.length; i += 2) args.set(Bun.argv[i], Bun.argv[i + 1]);

const repeats = Number(args.get("--repeats") ?? 1);
const exhaustiveness = Number(args.get("--exhaustiveness") ?? 8);
const cpuList = (args.get("--cpus") ?? "1,2,4,8").split(",").map(Number);
const seed = Number(args.get("--seed") ?? 20260717);
const timestamp = new Date().toISOString().replaceAll(":", "-").replace(".", "-");
const runDir = join(root, "results", timestamp);
await mkdir(runDir, { recursive: true });

const runs: Run[] = [];
for (const cpu of cpuList) {
  for (let repeat = 1; repeat <= repeats; repeat++) {
    const output = join(runDir, `pose-cpu${cpu}-run${repeat}.pdbqt`);
    const log = join(runDir, `vina-cpu${cpu}-run${repeat}.log`);
    const command = [
      join(root, "bin", "vina"),
      "--receptor", join(root, "data", "1iep_receptor.pdbqt"),
      "--ligand", join(root, "data", "1iep_ligand.pdbqt"),
      "--config", join(root, "data", "1iep_receptor.box.txt"),
      "--cpu", String(cpu),
      "--exhaustiveness", String(exhaustiveness),
      "--seed", String(seed),
      "--out", output,
    ];

    console.log(`Running CPU=${cpu}, repeat=${repeat}/${repeats}...`);
    const started = performance.now();
    const process = Bun.spawn(command, { stdout: "pipe", stderr: "pipe" });
    const [stdout, stderr, exitCode] = await Promise.all([
      new Response(process.stdout).text(),
      new Response(process.stderr).text(),
      process.exited,
    ]);
    const seconds = (performance.now() - started) / 1000;
    await Bun.write(log, stdout + stderr);
    const match = stdout.match(/^\s*1\s+(-?\d+(?:\.\d+)?)/m);
    const affinityKcalMol = match ? Number(match[1]) : null;
    runs.push({ cpu, repeat, seconds, affinityKcalMol, exitCode });
    console.log(`  ${seconds.toFixed(3)} s, affinity ${affinityKcalMol ?? "not parsed"} kcal/mol`);
    if (exitCode !== 0) throw new Error(`Vina failed; see ${log}`);
  }
}

const baselines = runs.filter((run) => run.cpu === cpuList[0]);
const baselineMedian = median(baselines.map((run) => run.seconds));
const summary = cpuList.map((cpu) => {
  const matching = runs.filter((run) => run.cpu === cpu);
  const medianSeconds = median(matching.map((run) => run.seconds));
  return {
    cpu,
    medianSeconds,
    speedup: baselineMedian / medianSeconds,
    efficiencyPercent: (baselineMedian / medianSeconds / cpu) * 100,
    bestAffinityKcalMol: Math.min(...matching.map((run) => run.affinityKcalMol ?? Infinity)),
  };
});

function median(values: number[]) {
  const sorted = [...values].sort((a, b) => a - b);
  const middle = Math.floor(sorted.length / 2);
  return sorted.length % 2 ? sorted[middle] : (sorted[middle - 1] + sorted[middle]) / 2;
}

const metadata = {
  generatedAt: new Date().toISOString(),
  vinaVersion: "1.2.7",
  system: "Apple M3 MacBook Air, 8 cores (4 performance + 4 efficiency), 8 GB RAM",
  receptor: "PDB 1IEP c-Abl kinase",
  ligand: "imatinib",
  exhaustiveness,
  seed,
  repeats,
  cpuList,
  runs,
  summary,
};
await Bun.write(join(runDir, "benchmark.json"), JSON.stringify(metadata, null, 2) + "\n");

const csv = [
  "cpu,median_seconds,speedup,efficiency_percent,best_affinity_kcal_mol",
  ...summary.map((row) => [
    row.cpu,
    row.medianSeconds.toFixed(3),
    row.speedup.toFixed(3),
    row.efficiencyPercent.toFixed(1),
    row.bestAffinityKcalMol.toFixed(3),
  ].join(",")),
].join("\n") + "\n";
await Bun.write(join(runDir, "summary.csv"), csv);

console.log("\n" + csv);
console.log(`Results: ${runDir}`);
