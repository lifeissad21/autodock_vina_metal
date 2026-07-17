import { existsSync } from "node:fs";
import { mkdir, readFile, writeFile } from "node:fs/promises";
import { basename, join, resolve } from "node:path";

const repo = resolve(import.meta.dir, "..");
const root = join(repo, "benchmark-100");
const vina = join(repo, "bin/vina");
const metal = join(repo, "metal-prototype/.build/release/VinaMetal");
const resultsRoot = join(root, "results");
const manifest = JSON.parse(await readFile(join(root, "manifest.json"), "utf8"));
const force = Bun.argv.includes("--force");
const limitAt = Bun.argv.indexOf("--limit");
const limit = limitAt >= 0 ? Number(Bun.argv[limitAt + 1]) : manifest.cases.length;
const seedBase = 20260717;
const exhaustiveness = 8;
const cpuThreads = 8;
const tty = process.stdout.isTTY;
const color = (code: number, text: string) => tty ? `\x1b[${code}m${text}\x1b[0m` : text;
const cyan = (text: string) => color(36, text), green = (text: string) => color(32, text);
const yellow = (text: string) => color(33, text), dim = (text: string) => color(2, text);
const benchmarkStarted = performance.now();

function duration(seconds: number) {
  if (seconds < 60) return `${seconds.toFixed(1)}s`;
  return `${Math.floor(seconds / 60)}m ${Math.round(seconds % 60)}s`;
}

function progress(position: number, total: number, completed = Math.floor(position)) {
  const width = 30;
  const exact = Math.max(0, Math.min(width, width * position / total));
  const filled = Math.floor(exact), partial = Math.floor((exact - filled) * 8);
  const fractions = ["", "▏", "▎", "▍", "▌", "▋", "▊", "▉"];
  const empty = width - filled - (partial ? 1 : 0);
  return `${cyan("█".repeat(filled) + fractions[partial])}${dim("░".repeat(empty))} ${String(completed).padStart(String(total).length)}/${total}`;
}

type Result = {
  id: string; target: string; ligand: string; seed: number; heavyAtoms: number; torsions: number;
  cpuSeconds: number; cpuScore: number; metalProcessSeconds: number; metalSearchSeconds: number;
  vinaRefinementSeconds: number; hybridScore: number; scoreDelta: number; poseRmsd: number | null;
};

async function run(command: string, args: string[], animate?: (seconds: number, frame: number) => void) {
  const started = performance.now();
  const process = Bun.spawn([command, ...args], { cwd: repo, stdout: "pipe", stderr: "pipe" });
  let frame = 0;
  const timer = tty && animate ? setInterval(() => animate((performance.now() - started) / 1000, frame++), 80) : undefined;
  const [exitCode, stdout, stderr] = await Promise.all([
    process.exited,
    new Response(process.stdout).text(),
    new Response(process.stderr).text(),
  ]);
  if (timer) clearInterval(timer);
  const wallSeconds = (performance.now() - started) / 1000;
  if (exitCode !== 0) {
    const tail = (stdout + stderr).trim().split("\n").slice(-16).join("\n");
    throw new Error(`${basename(command)} failed (${exitCode ?? "signal"})\n${tail}`);
  }
  return { stdout, stderr, wallSeconds };
}

function parse(text: string, pattern: RegExp, label: string) {
  const match = text.match(pattern);
  if (!match) throw new Error(`Could not parse ${label}`);
  return Number(match[1]);
}

async function heavyCoordinates(path: string) {
  const result: [number, number, number][] = [];
  let inFirstModel = false, sawModel = false;
  for (const line of (await readFile(path, "utf8")).split(/\r?\n/)) {
    if (line.startsWith("MODEL")) { if (sawModel) break; sawModel = inFirstModel = true; continue; }
    if (inFirstModel && line.startsWith("ENDMDL")) break;
    if (!line.startsWith("ATOM") && !line.startsWith("HETATM")) continue;
    const type = line.trim().split(/\s+/).at(-1)!;
    if (type === "H" || type === "HD") continue;
    result.push([Number(line.slice(30, 38)), Number(line.slice(38, 46)), Number(line.slice(46, 54))]);
  }
  return result;
}

async function rmsd(a: string, b: string) {
  const [left, right] = await Promise.all([heavyCoordinates(a), heavyCoordinates(b)]);
  if (!left.length || left.length !== right.length) return null;
  const squared = left.reduce((sum, point, i) => sum + point.reduce((s, value, axis) => s + (value - right[i][axis]) ** 2, 0), 0);
  return Math.sqrt(squared / left.length);
}

await mkdir(resultsRoot, { recursive: true });
const selected = manifest.cases.slice(0, limit);
const results: Result[] = [];
const spinner = ["⠋","⠙","⠹","⠸","⠼","⠴","⠦","⠧","⠇","⠏"];
function animateCase(index: number, id: string, phase: "CPU" | "Metal", seconds: number, frame: number) {
  const phaseStart = phase === "CPU" ? 0 : 0.48;
  const phaseSpan = phase === "CPU" ? 0.46 : 0.50;
  const smooth = 1 - Math.exp(-seconds / (phase === "CPU" ? 3.5 : 2.5));
  const position = index + phaseStart + phaseSpan * smooth;
  const label = phase === "CPU" ? "stock Vina" : "Metal + Vina refinement";
  process.stdout.write(`\x1b[2K\r${progress(position, selected.length, index)}  ${yellow(spinner[frame % spinner.length])} ${id.padEnd(10)} ${label}  ${dim(duration(seconds))}`);
}
console.log(`\n${cyan("VINA METAL BENCHMARK")}  ${dim(manifest.benchmark_version)}`);
console.log(`${dim("Cases")} ${selected.length}  ${dim("Proteins")} ${manifest.targets.length}  ${dim("Exhaustiveness")} ${exhaustiveness}  ${dim("CPU threads")} ${cpuThreads}`);
console.log(`${dim("Results")} ${resultsRoot}\n`);

for (const [index, test] of selected.entries()) {
  const id = `${test.target}-${String(index % manifest.ligands_per_target + 1).padStart(2, "0")}`;
  const caseDir = join(resultsRoot, id);
  const resultPath = join(caseDir, "result.json");
  if (!force && existsSync(resultPath)) {
    results.push(JSON.parse(await readFile(resultPath, "utf8")));
    console.log(`${progress(index + 1, selected.length)}  ${dim("↷")} ${id.padEnd(10)} ${dim("resumed")}`);
    continue;
  }
  await mkdir(caseDir, { recursive: true });
  const seed = seedBase + index;
  const cpuOut = join(caseDir, "vina-cpu.pdbqt");
  if (tty) process.stdout.write(`${progress(index, selected.length)}  ${yellow("●")} ${id.padEnd(10)} stock Vina...\r`);
  const cpu = await run(vina, ["--receptor", test.receptor, "--ligand", test.pdbqt, "--config", test.config,
    "--cpu", String(cpuThreads), "--exhaustiveness", String(exhaustiveness), "--seed", String(seed),
    "--num_modes", "1", "--out", cpuOut], (seconds, frame) => animateCase(index, id, "CPU", seconds, frame));
  await writeFile(join(caseDir, "vina-cpu.log"), cpu.stdout + cpu.stderr);

  const metalRaw = join(caseDir, "metal-raw.pdbqt"), metalFinal = join(caseDir, "metal-vina-final.pdbqt");
  if (tty) process.stdout.write(`${progress(index, selected.length)}  ${yellow("●")} ${id.padEnd(10)} Metal search + Vina refinement...\r`);
  const hybrid = await run(metal, ["--dock", "--flexible", "--maps", test.maps, "--ligand", test.pdbqt,
    "--output", metalRaw, "--exhaustiveness", String(exhaustiveness), "--num-modes", "1", "--seed", String(seed),
    "--vina-binary", vina, "--vina-receptor", test.receptor, "--vina-config", test.config, "--vina-output", metalFinal],
    (seconds, frame) => animateCase(index, id, "Metal", seconds, frame));
  await writeFile(join(caseDir, "metal-hybrid.log"), hybrid.stdout + hybrid.stderr);

  const cpuScore = parse(cpu.stdout, /^\s*1\s+(-?\d+(?:\.\d+)?)/m, "stock Vina score");
  const hybridScore = parse(hybrid.stdout, /Authoritative Vina 1\.2\.7 score:\s*(-?[\d.]+)/, "hybrid score");
  const result: Result = {
    id, target: test.target, ligand: test.ligand, seed, heavyAtoms: test.heavy_atoms, torsions: test.branches,
    cpuSeconds: cpu.wallSeconds, cpuScore, metalProcessSeconds: hybrid.wallSeconds,
    metalSearchSeconds: parse(hybrid.stdout, /Metal flexible docking:.*?,\s*([\d.]+) ms/, "Metal search time") / 1000,
    vinaRefinementSeconds: parse(hybrid.stdout, /Vina finalization:\s*([\d.]+) ms/, "Vina refinement time") / 1000,
    hybridScore, scoreDelta: hybridScore - cpuScore, poseRmsd: await rmsd(cpuOut, metalFinal),
  };
  await writeFile(resultPath, JSON.stringify(result, null, 2) + "\n");
  results.push(result);
  const elapsed = (performance.now() - benchmarkStarted) / 1000;
  const eta = elapsed / (index + 1) * (selected.length - index - 1);
  if (tty) process.stdout.write("\x1b[2K\r");
  console.log(`${progress(index + 1, selected.length)}  ${green("✓")} ${id.padEnd(10)} CPU ${duration(result.cpuSeconds).padStart(7)}  Hybrid ${duration(result.metalProcessSeconds).padStart(7)}  ${green((result.cpuSeconds/result.metalProcessSeconds).toFixed(2) + "x")}  Δ ${result.scoreDelta >= 0 ? "+" : ""}${result.scoreDelta.toFixed(3)}  ${dim(`ETA ${duration(eta)}`)}`);
}

const group = (rows: Result[]) => ({
  cases: rows.length,
  cpuSeconds: rows.reduce((s, r) => s + r.cpuSeconds, 0),
  hybridSeconds: rows.reduce((s, r) => s + r.metalProcessSeconds, 0),
  metalSearchSeconds: rows.reduce((s, r) => s + r.metalSearchSeconds, 0),
  meanAbsoluteScoreDelta: rows.length ? rows.reduce((s, r) => s + Math.abs(r.scoreDelta), 0) / rows.length : 0,
  maxAbsoluteScoreDelta: rows.length ? Math.max(...rows.map(r => Math.abs(r.scoreDelta))) : 0,
});
const targets = Object.fromEntries(manifest.targets.map((target: string) => [target, group(results.filter(r => r.target === target))]));
const totals = group(results);
const completedAt = new Date().toISOString();
const summary = { benchmarkVersion: manifest.benchmark_version, completedAt, seedBase, exhaustiveness, cpuThreads, totals, targets, results };
const header = "id,target,ligand,heavy_atoms,torsions,seed,cpu_seconds,hybrid_seconds,metal_search_seconds,vina_refinement_seconds,speedup,cpu_score,hybrid_score,score_delta,pose_rmsd";
const csv = [header, ...results.map(r => [r.id,r.target,r.ligand,r.heavyAtoms,r.torsions,r.seed,r.cpuSeconds,r.metalProcessSeconds,r.metalSearchSeconds,r.vinaRefinementSeconds,r.cpuSeconds/r.metalProcessSeconds,r.cpuScore,r.hybridScore,r.scoreDelta,r.poseRmsd ?? ""].join(","))].join("\n") + "\n";
const rows = Object.entries(targets).map(([target, value]: any) => `| ${target.toUpperCase()} | ${value.cases} | ${value.cpuSeconds.toFixed(2)} | ${value.hybridSeconds.toFixed(2)} | ${(value.cpuSeconds/value.hybridSeconds).toFixed(2)}x | ${value.meanAbsoluteScoreDelta.toFixed(3)} | ${value.maxAbsoluteScoreDelta.toFixed(3)} |`).join("\n");
const report = `# ${manifest.benchmark_version} results\n\nCompleted ${completedAt} on Apple M3. Stock AutoDock Vina 1.2.7 used ${cpuThreads} CPU threads; the hybrid used Metal search followed by authoritative Vina 1.2.7 explicit-receptor refinement. Both used exhaustiveness ${exhaustiveness}, one output mode, and matching per-case seeds. Map preparation is excluded.\n\n| Protein | Dockings | CPU Vina s | Metal + Vina s | Speedup | Mean abs score delta | Max abs score delta |\n|---|---:|---:|---:|---:|---:|---:|\n${rows}\n\n**Total:** ${results.length} dockings; CPU Vina ${totals.cpuSeconds.toFixed(2)} s; Metal + Vina ${totals.hybridSeconds.toFixed(2)} s; ${(totals.cpuSeconds/totals.hybridSeconds).toFixed(2)}x end-to-end speedup. Metal kernels totaled ${totals.metalSearchSeconds.toFixed(2)} s. Mean absolute final-score delta was ${totals.meanAbsoluteScoreDelta.toFixed(3)} kcal/mol; maximum was ${totals.maxAbsoluteScoreDelta.toFixed(3)} kcal/mol.\n\nScore delta is hybrid minus stock Vina. Both final scores are calculated by official Vina 1.2.7, but independent stochastic searches can select different minima. Pose RMSD in the CSV is direct atom-order heavy-atom RMSD and is not symmetry corrected.\n`;
await Promise.all([
  writeFile(join(resultsRoot, "summary.json"), JSON.stringify(summary, null, 2) + "\n"),
  writeFile(join(resultsRoot, "results.csv"), csv),
  writeFile(join(resultsRoot, "REPORT.md"), report),
]);
console.log(`\n${green("Benchmark complete")} in ${duration((performance.now() - benchmarkStarted) / 1000)}`);
console.log(`CPU Vina ${duration(totals.cpuSeconds)}  •  Metal + Vina ${duration(totals.hybridSeconds)}  •  ${green((totals.cpuSeconds/totals.hybridSeconds).toFixed(2) + "x speedup")}`);
console.log(`Mean |Δ score| ${totals.meanAbsoluteScoreDelta.toFixed(3)} kcal/mol  •  Report ${join(resultsRoot, "REPORT.md")}\n`);
