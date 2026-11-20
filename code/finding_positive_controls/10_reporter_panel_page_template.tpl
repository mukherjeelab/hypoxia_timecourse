<title>eIF3d Reporter Panel</title>
<meta name="claude:viewport-fit-sides" content="cover">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:wght@400;500;600&family=IBM+Plex+Mono:wght@400;500&family=Space+Grotesk:wght@500;600&display=swap">
<style>
/* Layout: one reading column; each table scrolls on its own; gene column sticky */
:root {
  --bg: #fbfbf9; --surface: #ffffff; --fg: #1d2125; --muted: #6a717a; --rule: #e2e4e1;
  --target: #E69F00; --pos: #0072B2; --neg: #1d2125;
  --sig: #1d2125; --nonsig: #9aa0a6;
  --excl: #e6f1ea; --excl-fg: #23663f; --buff: #fdf1df; --buff-fg: #8a5a00;
  --fwd: #e8eef6; --fwd-fg: #2f5585; --ns: #eff0ee; --ns-fg: #5d646c;
  --display: "Space Grotesk", system-ui, sans-serif;
  --body: "IBM Plex Sans", system-ui, sans-serif;
  --mono: "IBM Plex Mono", ui-monospace, Menlo, monospace;
}
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg: #15181b; --surface: #1c2024; --fg: #e7e9ea; --muted: #98a0a8; --rule: #30353a;
  --neg: #e7e9ea; --sig: #e7e9ea; --nonsig: #6f767d;
  --excl: #1f3328; --excl-fg: #8fd1a8; --buff: #3a2e17; --buff-fg: #f0c27a;
  --fwd: #1f2a3a; --fwd-fg: #9dbbe3; --ns: #262a2e; --ns-fg: #a3aab1; color-scheme: dark } }
:root[data-theme="dark"] {
  --bg: #15181b; --surface: #1c2024; --fg: #e7e9ea; --muted: #98a0a8; --rule: #30353a;
  --neg: #e7e9ea; --sig: #e7e9ea; --nonsig: #6f767d;
  --excl: #1f3328; --excl-fg: #8fd1a8; --buff: #3a2e17; --buff-fg: #f0c27a;
  --fwd: #1f2a3a; --fwd-fg: #9dbbe3; --ns: #262a2e; --ns-fg: #a3aab1; color-scheme: dark }
body { background: var(--bg); color: var(--fg); font-family: var(--body); font-size: 15px; line-height: 1.5; }
.wrap { max-width: 1120px; margin: 0 auto; padding-block: 40px 64px;
  padding-left: max(16px, env(safe-area-inset-left, 0px)); padding-right: max(16px, env(safe-area-inset-right, 0px));
  display: flex; flex-direction: column; gap: 40px; }
header { display: flex; flex-direction: column; gap: 10px; max-width: 70ch; }
h1 { font-family: var(--display); font-weight: 600; font-size: 2rem; margin: 0; text-wrap: balance; letter-spacing: -0.01em; }
h2 { font-family: var(--display); font-weight: 600; font-size: 1.25rem; margin: 0; text-wrap: balance; }
.lede { margin: 0; color: var(--muted); }
.legend { display: flex; flex-wrap: wrap; gap: 8px 20px; font-size: 0.85rem; color: var(--muted); }
.legend span { display: inline-flex; align-items: center; gap: 6px; }
section { display: flex; flex-direction: column; gap: 12px; }
section > p { margin: 0; max-width: 75ch; color: var(--muted); font-size: 0.9rem; }
.scroll { overflow-x: auto; border: 1px solid var(--rule); border-radius: 6px; background: var(--surface); min-width: 0; }
table { border-collapse: collapse; width: 100%; font-variant-numeric: tabular-nums; }
th, td { padding: 9px 12px; text-align: left; vertical-align: top; border-bottom: 1px solid var(--rule); white-space: nowrap; }
tbody tr:last-child th, tbody tr:last-child td { border-bottom: 0; }
thead th { font-size: 0.72rem; text-transform: uppercase; letter-spacing: 0.06em; color: var(--muted); font-weight: 500; background: var(--surface); }
thead tr.group th { font-family: var(--display); text-transform: none; letter-spacing: 0; font-size: 0.85rem; color: var(--fg); font-weight: 600; border-bottom: 0; padding-bottom: 0; }
th[scope="row"] { position: sticky; left: 0; background: var(--surface); font-weight: 600; z-index: 1; }
th[scope="row"] .role { display: block; font-size: 0.72rem; font-weight: 400; color: var(--muted); padding-left: 16px; }
.split { border-left: 1px solid var(--rule); }
.dot { display: inline-block; width: 9px; height: 9px; border-radius: 50%; margin-right: 7px; vertical-align: 1px; }
.dot.target { background: var(--target); } .dot.pos { background: var(--pos); } .dot.neg { background: var(--neg); }
.v { font-family: var(--mono); font-size: 0.88rem; color: var(--nonsig); }
.v.sig { color: var(--sig); font-weight: 500; }
.p { font-family: var(--mono); font-size: 0.72rem; color: var(--muted); margin-left: 6px; }
.cat { font-size: 0.75rem; font-weight: 500; padding: 2px 8px; border-radius: 999px; }
.cat-exclusive { background: var(--excl); color: var(--excl-fg); } .cat-buffered { background: var(--buff); color: var(--buff-fg); }
.cat-forwarded { background: var(--fwd); color: var(--fwd-fg); } .cat-not { background: var(--ns); color: var(--ns-fg); }
.cat-te_only { background: var(--excl); color: var(--excl-fg); }
.tgt { white-space: normal; min-width: 14ch; font-size: 0.85rem; }
.cat-not_sig { background: var(--ns); color: var(--ns-fg); }
.muted { color: var(--muted); }
.part { display: flex; flex-direction: column; gap: 4px; padding-top: 20px; border-top: 2px solid var(--fg); }
.part p { margin: 0; color: var(--muted); font-size: 0.9rem; }
.eyebrow { font-size: 0.72rem; text-transform: uppercase; letter-spacing: 0.08em; color: var(--muted); font-weight: 500; }
.part-title { font-size: 1.5rem; }
.tgt .call { display: block; font-weight: 600; }
.tgt .rna { display: block; font-size: 0.78rem; color: var(--muted); margin-top: 2px; }
.notes-lede { margin: 0; color: var(--muted); font-size: 0.9rem; }
.notes { display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 16px 28px; }
.note { min-width: 0; } .note h3 { font-family: var(--display); font-size: 1rem; margin: 0 0 4px; }
.note p { margin: 0; font-size: 0.9rem; }
dl { display: grid; grid-template-columns: max-content 1fr; gap: 6px 16px; margin: 0; font-size: 0.88rem; max-width: 90ch; }
dt { font-weight: 600; } dd { margin: 0; color: var(--muted); }
@media (max-width: 560px) { dl { grid-template-columns: 1fr; } dd { margin-bottom: 8px; } h1 { font-size: 1.6rem; } }
</style>

<div class="wrap">
<header>
  <h1>eIF3d Reporter Panel</h1>
  <p class="lede">Nine candidate 5'UTR reporter genes: five eIF3d targets selected from hypoxia ribosome profiling, two positive controls from Roiuk et al. 2024, and two negative controls. ENOX2 was added last, from a screen for genes whose eIF3d effect is purely translational (Exclusive) in both lines and absent in normoxia. The first part is our own ribosome profiling of eIF3d (si3d) and eIF3e (si3e) knockdown in MDA-MB-231 and MCF7-SIX1. The second part is published data from other labs, kept separate: eIF3d and eIF3e knockdown in HeLa, and eIF3 binding.</p>
  <div class="legend">
    <span><span class="dot target"></span>Selected eIF3d target</span>
    <span><span class="dot pos"></span>Positive control (Roiuk et al.)</span>
    <span><span class="dot neg"></span>Negative control</span>
  </div>
</header>

<div class="part"><span class="eyebrow">Part 1</span><h2 class="part-title">Our data: MDA-MB-231 and MCF7-SIX1</h2></div>

<section>
  <h2>eIF3d knockdown, hypoxia 1 hr</h2>
  <p>log2 fold change, sictrl / si3d, with adjusted p. Positive means eIF3d knockdown lowers it, so positive TE means eIF3d promotes translation. Dark values are significant (padj &lt; 0.05); grey values are not.</p>
  <div class="scroll"><table>
    <thead>
      <tr class="group"><th></th><th colspan="4">MDA-MB-231</th><th colspan="4" class="split">MCF7-SIX1</th></tr>
      <tr><th>Gene</th><th>TE</th><th>Ribo</th><th>mRNA</th><th>Category</th><th class="split">TE</th><th>Ribo</th><th>mRNA</th><th>Category</th></tr>
    </thead>
    <tbody>
%%T1%%
    </tbody>
  </table></div>
</section>

<section>
  <h2>eIF3e knockdown, hypoxia 1 hr</h2>
  <p>Same layout, log2 fold change sictrl / si3e. An eIF3e target has TE log2 fold change above 0.5 with padj &lt; 0.05, the same rule used for eIF3d targets.</p>
  <div class="scroll"><table>
    <thead>
      <tr class="group"><th></th><th colspan="4">MDA-MB-231</th><th colspan="4" class="split">MCF7-SIX1</th><th class="split"></th></tr>
      <tr><th>Gene</th><th>TE</th><th>Ribo</th><th>mRNA</th><th>Category</th><th class="split">TE</th><th>Ribo</th><th>mRNA</th><th>Category</th><th class="split">eIF3e target? (TE)</th></tr>
    </thead>
    <tbody>
%%T4%%
    </tbody>
  </table></div>
</section>

<section>
  <h2>eIF3d knockdown, normoxia</h2>
  <p>TE log2 fold change, sictrl / si3d (padj). Compare with hypoxia above to see whether the eIF3d effect is specific to hypoxia.</p>
  <div class="scroll"><table>
    <thead><tr><th>Gene</th><th>MDA-MB-231 (1 + 4 hr, time-regressed)</th><th>MCF7-SIX1 (1 hr)</th></tr></thead>
    <tbody>
%%T2%%
    </tbody>
  </table></div>
</section>

<div class="part"><span class="eyebrow">Part 2</span><h2 class="part-title">Published data from other labs</h2><p>Not our experiments. Shown for independent support only.</p></div>

<section>
  <h2>HeLa eIF3d and eIF3e knockdown</h2>
  <p>Herrmannová et al., eLife 2024. Ribosome profiling in HeLa at steady state (normoxia). TE log2 fold change, sictrl / knockdown (padj), oriented the same way as our data.</p>
  <div class="scroll"><table>
    <thead><tr><th>Gene</th><th>si3d</th><th>si3e</th></tr></thead>
    <tbody>
%%T5%%
    </tbody>
  </table></div>
</section>

<section>
  <h2>eIF3 binding</h2>
  <p>"Not called" means the gene is absent from that study's hit list, not that binding was measured and ruled out. Detection follows expression, and none of these studies used our cell lines.</p>
  <div class="scroll"><table>
    <thead><tr><th>Gene</th><th>eIF3d cap binding<br>(subunit-seq, HEK293T)</th><th>Lee 2015 PAR-CLIP<br>(HEK293T)</th><th>De Silva 2021 PAR-CLIP<br>eIF3d band (Jurkat)</th><th>De Silva 2021<br>other eIF3 bands</th><th>Cate lab 2023 irCLIP<br>(NPCs)</th></tr></thead>
    <tbody>
%%T3%%
    </tbody>
  </table></div>
</section>

<section>
  <h2>Notes</h2>
  <p class="notes-lede">Mostly from our data; anything from a published study is labelled as such.</p>
  <div class="notes">
%%NOTES%%
  </div>
</section>

<section>
  <h2>Definitions and methods</h2>
  <dl>
    <dt>TE</dt><dd>Translation efficiency, ribosome footprints relative to mRNA, from a DESeq2 interaction model.</dd>
    <dt>Exclusive</dt><dd>Footprints and TE change; mRNA does not. A purely translational effect.</dd>
    <dt>Buffered</dt><dd>mRNA and TE both change (here mRNA rises on knockdown), so the protein-level effect is likely smaller than the TE value suggests.</dd>
    <dt>Forwarded</dt><dd>mRNA and footprints change together; TE does not.</dd>
    <dt>Not sig</dt><dd>No significant TE change. Categories use padj &lt; 0.05 with no fold-change cutoff.</dd>
    <dt>MDA-MB-231</dt><dd>Spike-in normalized, apeglm-shrunken fold changes. Knockdown lowers total mRNA per cell, so the typical gene sits near +0.4 for mRNA; an mRNA value near that is unremarkable.</dd>
    <dt>MCF7-SIX1</dt><dd>Not spike-in normalized. TE is apeglm-shrunken; ribo and mRNA are unshrunk fold changes from the same samples.</dd>
    <dt>eIF3e mRNA call</dt><dd>"Higher" or "lower" after knockdown means padj &lt; 0.05 and at least 0.25 log2 away from the median gene in that experiment. Measured against the median because MDA-MB-231 mRNA is spike-in normalized and the typical gene sits near +0.4.</dd>
    <dt>TE values are per gene</dt><dd>They reflect the mix of isoforms each cell line expresses. A reporter carries one 5'UTR, so check that it matches the expressed transcript (see notes).</dd>
  </dl>
</section>
</div>
