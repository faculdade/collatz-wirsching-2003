# E-010: central-cost zeros and the quantitative half of Wirsching's (?3)

Related hypotheses: H-003 (dated section of 2026-08-09), H-012, H-013,
H-014 (opened from this experiment's quantile run), H-015 (this
experiment is its entire evidence base) and H-019 (opened from the
ell=19/20 break below).

Wirsching (2003) reduces uniform positive predecessor density to a chain
of five conditions. Conjecture 2 is the step `(?4) => (?3)`, where

```text
(?3)  g_ell(k_ell, a) >= mu * gbar_ell(k_ell)
      for every unit a, every ell >= ell_0, and every sequence
      (k_ell) with |ell - k_ell| <= delta*sqrt(ell),
```

and `(?4)` is a statement about the one-dimensional averaging operator
`W_3` that never mentions the generators. Since `k_ell = ell` is an
admissible sequence for every `delta > 0`, a single unit residue with
`g_ell(ell,a) = 0` at one level `ell >= ell_0` falsifies `(?3)`.

This experiment separates the two independent halves of `(?3)`.

## `central_zeros.py`: the support half

Exact Boolean support of `g_ell(.,a)` for every residue modulo `3^ell`
and every cost up to a ceiling, packed one bit per cost into one word
per residue. Reports the zero set at cost exactly `ell`, the least cost
whose support covers every unit, and the coherent subtree of residues
whose truncations are central-cost zeros at every earlier level.

```sh
python3 central_zeros.py --max-ell 18
```

Through `ell = 18` the zero set is never empty and the least fully
covering cost is `ell + 5` for every `ell` from 10 on, extending the
pattern E-004 saw through 16. At `ell = 18` there are 11,540,739 zeros
among 258,280,326 units and 734,754 coherent ones. A coherent subtree
that stayed nonempty at every level would produce, by the inverse limit
of nonempty finite sets, a 3-adic `alpha` that is a central-cost zero at
every level, refuting `(?3)` and `(?2)` at that `alpha`. Its growth
ratio falls from 2.17 to 1.79 over `ell = 15..18`, which by an
endpoint-slope read of that window (-0.1280/level) put the crossing of
1 near `ell = 24`. Updated 2026-08-16 with the two new levels: the
ratios over `ell = 17..20` are 1.94, 1.79, 1.633, 1.553. Applying the
SAME method to that window gives -0.1287/level and a crossing at
`ell = 24.3`, i.e. essentially unchanged, and 4-point OLS gives 24.0.
(A previous version of this paragraph reported "roughly `ell = 27`"
and said the decline had slowed; both were wrong. The 27 came from
reading only the last single step, 19 to 20, under a different method
from the one that produced 24, and the window-average decline is if
anything marginally faster, not slower. Corrected in critique round 3.)
This is extrapolation either way, and not a claim about what actually
happens.

**Update, `ell = 19` and `ell = 20` (2026-08-15):** the `ell + 5`
pattern above does not continue. `first_full - ell` reads `+4` at
`ell = 19` and `+7` at `ell = 20`, the first break in that pattern
since it started holding steady at `ell = 10`. `frac_zero` is
non-monotone across the same two levels (0.0447 at 18, down to 0.0315
at 19, back up to 0.0381 at 20), the first break since its long,
otherwise unbroken decline set in (it was also non-monotone once
before that, rising from 0.5000 at `ell = 1` to 0.6667 at `ell = 2`,
an opening transient at the very start of the series, not part of the
established decreasing regime this update is comparing against). The
coherent-subtree
count keeps growing (1,199,856 at 19, 1,863,109 at 20) but its growth
ratio keeps falling (1.633, then 1.553, continuing the decline from
1.79 at `ell = 18`), consistent with, not proof of, the same eventual-
extinction extrapolation. Reported as observed; no claim about what
either irregularity means beyond these two data points (Rule 11).

The packed table is audited against the independent backward predicate
of E-007 at every level through 7, and the smallest coherent witnesses
of each level are rechecked against the same predicate at every earlier
level.

## `central_ratio.py`: the quantitative half

Exact integer counts from the same recursion, reporting

```text
min_{a in S} g_ell(ell+d, a) / gbar_ell(ell+d)
```

for a range of offsets `d`, and separately for a fixed set of integers
`a` that does not grow with `ell`. Every row total is checked against
the independent count of bounded compositions of `k` with capacities
`2, 6, 18, ..., 2*3^(ell-1)`, which is Wirsching's identity
`2*3^(ell-1)*gbar_ell(k)`.

```sh
python3 central_ratio.py --max-ell 15
python3 central_ratio.py --max-ell 16 --offsets 0 5 12 --sqrt-multiples 1
python3 central_ratio.py --max-ell 16 --offsets 5 --sqrt-multiples 1 \
    --fixed-offset 5 --random-fixed 1458 --random-bound 2187
```

Minimum over all units, at offsets where the support is already
complete, decays geometrically in `ell`: at `d = +12` it falls from
0.4805 at `ell = 6` to 0.2555 at `ell = 16`, a factor near 0.94 per
level that is the same on both halves of the range. The offset would
have to grow linearly in `ell` to compensate, and Wirsching's window
allows only `d <= delta*sqrt(ell)`. So the infimum over `Z_3^x` that
`(?3)` demands stabilises nowhere in the window, not only at `k = ell`
for lack of support.

Minimum over a fixed set of integers behaves differently. At `d = +5`,
over all 486 units below `3^6`, the ratio reads 0.1019, 0.1472, 0.1341,
0.1377, 0.1306, 0.1397, 0.1038 for `ell = 10..16`: flat. Over all 1458
units below `3^7` it reads 0.102, 0.074, 0.107, 0.079, 0.116, 0.108,
0.096 for `ell = 10..16`, and over all 4374 units below `3^8` it reads
0.0509, 0.0736, 0.0805, 0.0787, 0.1016, 0.0967, 0.0878: also flat.
Over all 54 units below `3^4` it reads 0.3056, 0.2576, 0.1877, 0.2361,
0.1742, 0.1720, 0.1677 for the same levels. **This one does decline**:
five of its six transitions fall, and it loses 45% from `ell = 10` to
`ell = 16`. That is different from the other three fixed sets above,
which are genuinely flat, and it is the smallest set here (a minimum
over 54 values, so also the noisiest). Two earlier summaries of this
line were both wrong: "sits near 0.22" (the middle of the range, not a
value it settles at) and then "no trend" (indefensible against a
45% monotone-ish fall). Corrected 2026-08-16 after the second was
caught in critique round 3.

Whether the smallest fixed set declining while larger ones stay flat
means anything is not settled here, and H-013, which drew the
"fixed sets do not deteriorate" conclusion, is closed
(`closed-inconclusive`) and is not reopened by this note (Rule 8d).
The observation is recorded so it is not lost; a session wanting to
pursue it should open a hypothesis rather than edit H-013's closure.

Each of these sets is exhaustive, not sampled: there are exactly 54,
486, 1458 and 4374 non-multiples of 3 below `3^4`, `3^6`, `3^7` and
`3^8` respectively, so `--random-fixed N --random-bound 3^m` with
those matched pairs enumerates the whole set. Raw data for the 54/486/
4374 series in `ratio_fixedsets_ell16.log` (generated 2026-08-16
during critique, after these three series were found to be quoted with
no committed output backing them; every previously quoted digit
reproduced exactly). The 1458 series is the one carried by
`ratio_ell18_output.txt`.

The minimum over all units modulo `3^ell` is the minimum over integers
`a < 3^ell`, and that is not evidence about `liminf_ell` at any fixed
`a`, because the residue attaining it has an integer representative of
size `3^ell`. What the tables establish is that on the three larger
fixed sets (486, 1458, 4374 integers) the ratio does not deteriorate
with `ell`, while the all-unit minimum decays; that decay is the growth
of the index set. The smallest set (54 integers) is the exception noted
above and does fall. Separating that from a genuine decay of the
left tail needs quantiles rather than the minimum, which is the next
step recorded in H-013.

## `central_quantiles.py`: is the minimum's decay extreme-value statistics?

Fixed quantiles of the same distribution `g_ell(k,a)/gbar_ell(k)` over
unit residues, `k = ell + d`, run through `ell = 16`:

```sh
python3 central_quantiles.py --max-ell 16 --offsets 5 12 \
    --quantiles 1e-4 1e-3 1e-2 1e-1 0.5
```

A quantile does not carry the extreme-value confound of a minimum: the
value below which a fixed *fraction* of the population falls does not
mechanically shrink just because the population grows. H-013 states the
test directly: if a low quantile stays flat in `ell` while the minimum
falls, the minimum's decay is extreme-value statistics of a growing
index set; if the quantile also falls, the left tail is deteriorating
for real.

The two offsets return opposite verdicts by that test. At `d = +5`
(support complete from `ell = 10`), `q = 10^-3` reads 0.102, 0.110,
0.107, 0.118, 0.116, 0.107, 0.104 for `ell = 10..16`: flat, while the
minimum over the same rows is noisy and small (0.0197-0.0536, corrected
2026-08-16 from a stated "0.02-0.07" whose upper end exceeds every
value in the column) with no clear trend. The criterion reads extreme-value statistics. At `d = +12`
(support complete from `ell = 4`), `q = 10^-4` falls from 0.402 at
`ell = 10` to 0.276 at `ell = 16`, a factor near 0.94 per level, close
to the minimum's own decay rate at that offset. The criterion reads
real tail deterioration. Neither offset is more authoritative than the
other; both are legitimate members of Wirsching's window
`|ell - k_ell| <= delta*sqrt(ell)`, and the disagreement between them is
itself the finding, not noise to average away (H-014, opened from this
run).

Part of the `d = +12` decay is bulk drift rather than tail pinching:
the median (`q = 0.5`) also falls, from 1.011 at `ell = 10` to 0.953 at
`ell = 16`, so `gbar_ell` is pulled up by the right tail relative to
the bulk. Normalizing `q = 10^-4` by the median instead of by 1 still
decays at close to the same 0.94/level rate, so bulk drift is not the
whole story at `d = +12`. At `d = +5` the same normalization makes the
already-flat `q = 10^-4` column rise slightly (0.053 to 0.093 over
`ell = 10..16`), the opposite direction. The `q = 10^-4` column is the
4th order statistic at `ell = 10` (`n_units = 39366`, so `10^-4` picks
index 3) and is min-contaminated there; read its trend only from
`ell >~ 13`.

Full data (`ell`, `d`, `n_units`, quantiles `1e-4/1e-3/1e-2/1e-1/0.5`,
min) in `quantiles_ell16.log` in this folder.

## Resolving the `d`-sensitivity: an offset grid (H-014)

`central_quantiles.py --max-ell 16 --offsets 6 7 8 9 10 11` (plus
`--thresholds` and `--bucket-quantile`, added to the script alongside
this run) settles which offset is the outlier. Every `d` from 6 to 11
decays at a rate close to `d = +12`'s (0.93-0.96 per level on
`q = 10^-3`); `d = +5` is the one that stays flat. `d = +5` is not a
generic point in the window: it is exactly the first cost with
complete support at each level (`central_zeros.py`'s own boundary,
`ell + 5` for every `10 <= ell <= 18`), so its distribution is
dominated by residues whose count has only just left zero, a boundary
regime distinct from `d >= 6` where support has been complete for
several levels already. Full data in `quantiles_grid_ell16.log`.

## What the population tail decay is made of: bottom-bucket composition

`central_quantiles.py` also reports, for the bottom `10^-3` bucket of
the unit population at each level, what fraction of those residues has
an integer representative below `3^8 = 6561` (a "small", fixed-in-`ell`
integer by the standard of the fixed-set columns above). That fraction
falls by roughly three orders of magnitude at every offset tested
(`d = 5` through `12`): at `ell = 9` it spans 0.23 to 0.62 across the
eight offsets (0.54 at `d=5`, 0.62 at `d=6`, 0.23-0.38 elsewhere), and
at `ell = 16` it spans 0.00007 to 0.00045. (Corrected 2026-08-16: an
earlier version of this sentence stated "about 0.5 at `ell = 9`" and
"0.0002-0.0007 at `ell = 16`", neither of which brackets the actual
per-offset values.) The population's low tail
is, at the levels that matter, almost entirely new residues whose
integer representative only exists from that level on. This is the
mechanism reconciling the population decay above with the flatness of
the three larger fixed sets: the group grows and keeps injecting bad
new residues.

**It does not establish that no small fixed `a` deteriorates**
(corrected 2026-08-16, round 4; an earlier version said exactly that,
and round 3's narrowing elsewhere in this file was not propagated
here). The exhaustive 54-integer set at the same offset falls 45% over
`ell = 10..16`, and because that set is fixed and finite, a falling
minimum over it means some `a < 81` does get worse. The bucket
composition and the 54-set are measuring different things and only the
larger sets support the flatness reading. This bucket-composition
column is part of every `central_quantiles.py` run's own output (not a
separate command); full data for every offset tested is in
`quantiles_ell16.log` (`d=5,12`) and `quantiles_grid_ell16.log`
(`d=6..11`) (regenerated 2026-08-15 during critique, after a prior
reference to a separate `quantiles_thresholds_ell16.log` was found to
point at a file that never existed in this repository; the data itself
was never lost, it always lived in these two files).

## Dynamic memory guard for open-ended `central_zeros.py` runs

`central_zeros.py` now measures its own actual RSS growth after each
completed level (via `/proc/self/statm`), projects the next level's
need as that delta times `--growth-factor` (default 3.3, padded above
the theoretical 3x since the modulus itself triples), and checks it
against currently available system RAM. If continuing would leave
less than `--reserve-gib` (default 20.0) GiB free system-wide, it
prints a `[stop]` message and exits the loop cleanly, keeping whatever
levels it already completed, rather than running until `MemoryError`
or degrading the whole machine. In RAM mode this lets a run be
launched with a generous `--max-ell` and trusted to stop itself.

**That guarantee does not hold in disk mode** (added later by
`a296a52`, and not described here until round 5 caught the omission).
`--scratch-dir` memmaps the packed word arrays only, and **disables the
RSS guard**, since resident memory no longer reflects the arrays; a
free-disk check replaces it, which does not bound the still-resident
Boolean arrays. At `ell = 21` those alone exceed this machine. See
`notes/H-019.md` for the numbers before launching anything past
`ell = 20`.

The same commit also changed the reporting of `first_full`: a bare
`none` used to mean both "no covering cost exists" and "the search ran
out of ceiling", and now `none` is reserved for the definitive case
while `>+N` marks the inconclusive one, decided by comparing against
the largest attainable cost `3^ell - 1 - ell`.

## Evidence (Rule 9a): `central_zeros.py`, extended to `ell = 20`

```
Command:      python3 central_zeros.py --max-ell 10
Commit:       43ce2f5
Date:         2026-08-17
Environment:  Linux, Python 3.12.3, numpy 2.5.1
Exit:         0
Output:       experiments/E-010-wirsching-conj2-central-zeros/zeros_regression_ell10.log
Checked:      producer. The regression run this README described as
              "verified correctness-preserving" had no committed output
              until now; critique round 5 raised it. Reproduces every
              level through ell=10, with first_full - ell equal to +6 at
              ell=8,9 and +5 at ell=10.

Command:      python3 central_ratio.py --max-ell 8
              python3 central_ratio.py --max-ell 8 --scratch-dir <tmp>
Commit:       43ce2f5
Date:         2026-08-17
Environment:  Linux, Python 3.12.3, numpy 2.5.1
Exit:         0 (both)
Output:       experiments/E-010-wirsching-conj2-central-zeros/ratio_ram_vs_disk_ell8.log
Checked:      producer, by diffing the two outputs after removing the
              scratch-directory preflight lines that exist only in disk
              mode. They agree line for line. The README asserted this
              byte-for-byte agreement with no retained comparison until
              critique round 5 raised it.

Command:      python3 weight_asymmetry.py --ell 4 --k 5
Commit:       8911937
Date:         2026-08-17
Environment:  Linux, Python 3.12.3
Exit:         0
Output:       experiments/E-010-wirsching-conj2-central-zeros/weight_asymmetry_ell4_k5.log
Checked:      producer, against main.tex's Remark 5.13. S_4(5,a) = 5/4 and
              7/4 on the two unit classes, so their sum is 3 and the Haar
              mean is 3/2, both asserted by the script. Added because
              critique round 3 found this pair quoted as an exact
              evaluation with no run behind it; the committed logs covered
              only S_3(3,.) and S_4(4,.).

Command:      python3 central_zeros.py --max-ell 30
Commit:       cfc63a1d21ebbfaca9bd93c6787d0715984d2975 (the tree the run
              actually executed against: the output file's mtime and its
              own per-level timers put the start at or before 11:47:38,
              while ab2b13b was authored 11:49:31. Corrected in critique
              round 3 from ab2b13b, the recording commit. central_zeros.py
              is byte-identical at both, so reproducibility is unaffected)
Date:         2026-08-15
Environment:  Linux 7.0.0-28-generic x86_64, Python 3.12.3, 62 GiB RAM
Exit:         0
Output:       experiments/E-010-wirsching-conj2-central-zeros/zeros_extended_output.txt
Checked:      producer (main session); every level through ell=18 matches
              the previously documented numbers exactly (independent
              reproduction in this repository); ell=19/20 are new, not
              previously computed here; internal audit against the
              backward predicate ran through ell=7 as usual. Critiqued in
              H-015 rounds 2 and 3 (see notes/H-015.md); round 3's R3-03 corrected an extrapolation
              drawn from this file's coherent counts. Rounds 4 and 5
              also drew on this file.
```

Self-stopped cleanly at `ell = 20` via the dynamic RAM guard (`level 21
projected to need ~59.5 GiB more, but only 26.3 GiB is available`), not
an error. Reproduces the `ell = 18` figures cited above exactly
(11,540,739 zeros, 734,754 coherent) and adds two new levels; see the
"Update, ell = 19 and ell = 20" note above for what changed there.

Verified correctness-preserving on 2026-08-15: `--max-ell 10` produces
the same computation and passes every internal audit assertion
(against the independent backward predicate) as before the change;
the new per-row `rss=X.XXGiB` note is additive only.

## Evidence (Rule 9a): `central_quantiles.py` logs, regenerated

The three quantile log files this document cites (`quantiles_ell16.log`,
`quantiles_grid_ell16.log`, previously also a non-existent
`quantiles_thresholds_ell16.log`, see the note above) were referenced
but never committed to this repository. Regenerated 2026-08-15 during
critique; every number checked against this document's own prose
(the flat `d=+5` quantile band, the `d=+12` 0.402-to-0.276 decay, the
bucket-composition fractions) matches exactly.

```
Command:      python3 central_quantiles.py --max-ell 16 --offsets 5 12
              --quantiles 1e-4 1e-3 1e-2 1e-1 0.5
              > quantiles_ell16.log
              python3 central_quantiles.py --max-ell 16 --offsets 6 7 8 9 10 11
              > quantiles_grid_ell16.log
Commit:       653cc5db0b741547ff559d84dcebc04791f48c40
Date:         2026-08-15
Environment:  Linux 7.0.0-28-generic x86_64, Python 3.12.3
Exit:         0 (both)
Output:       experiments/E-010-wirsching-conj2-central-zeros/quantiles_ell16.log,
              experiments/E-010-wirsching-conj2-central-zeros/quantiles_grid_ell16.log
Checked:      producer (main session), against the specific numbers already
              quoted in this document's "central_quantiles.py" and
              "Resolving the d-sensitivity" and "bottom-bucket composition"
              sections. Critiqued in H-015 rounds 2 and 3 (see
              notes/H-015.md); round 2's R2-01 corrected the
              bucket-composition interval quoted from these logs.
```

## Evidence (Rule 9a): the exhaustive fixed-set series

```
Command:      python3 central_ratio.py --max-ell 16 --offsets 5
              --sqrt-multiples 1 --fixed-offset 5
              --random-fixed {54,486,4374} --random-bound {81,729,6561}
              (three runs, one per set, concatenated into one log)
Commit:       397607e29cf42231ff3aaeea0f4725db12408c0f
Date:         2026-08-16
Environment:  Linux 7.0.0-28-generic x86_64, Python 3.12.3, 62 GiB RAM
Exit:         0 (all three)
Output:       experiments/E-010-wirsching-conj2-central-zeros/ratio_fixedsets_ell16.log
Checked:      producer (main session): every digit previously quoted in
              this document for the 486 and 4374 series reproduced
              exactly; the 54 series' quoted summary was corrected
              against the actual values. Critiqued in H-015 round 3,
              which corrected the 54-series characterization drawn from
              this log (see notes/H-015.md).
```

## Extending the fixed-integer test to `ell = 17`

`central_ratio.py --max-ell 17 --offsets 5 12 --sqrt-multiples 1
--fixed-offset 5 --random-fixed 1458 --random-bound 2187` (peak
38.5 GiB by the corrected model, `k_max = 29`: 28.87 GiB for the
ell=17 table plus the still-live ell=16 one. An earlier version quoted
"~39 GiB checked live against free RAM"; that observation came from the
same run whose log was never committed, so it is not backed by anything
in this repository and has been replaced by the derived figure.
Corrected in critique round 3, the same defect class round 2 fixed for
that run's timing.) The
minimum over the same exhaustive 1458 integers below `3^7` at `d = +5`
reads 0.102, 0.074, 0.107, 0.079, 0.116, 0.108, 0.096, 0.089 for
`ell = 10..17`: eight levels, no bend at the new one either. The raw
log this originally cited (`ratio_ell17.log`) was never committed to
this repository (found during critique, 2026-08-15); the same ell=17
row, however, is exactly reproduced as a subset of the `ell = 1..18`
table in `ratio_ell18_output.txt` (the disk-backed rerun below), whose
`d=+5` fixed-integer column at `ell=17` reads the identical `0.0892`,
so this specific data is not lost, just backed by a differently-named
file than first stated.

## Escalation: does H-166's cascade transfer here? (Regra 11b)

A structurally different result in this project (H-166, unrelated to
Wirsching's `g_ell`) proves `min_u N_ell(u)/N_(ell-1)(...)` is a convex
combination of the previous level's ratios, giving monotonicity and a
bound certified at every higher level from one finite computation. An
external model (Codex, high reasoning effort) was asked whether the
same mechanism applies to `R_ell(k,a) = g_ell(k,a)/gbar_ell(k)`. It
does not, with a derivation and an explicit counterexample rather than
a bare assertion:

- Dividing the `extend` recursion by `gbar_ell(k)` gives weights whose
  sum over the admissible `j`'s for a fixed `a`, `S_ell(k,a)`, is **not**
  identically 1, and **not on Haar average either**: the identity is
  `S_ell(k, a=1) + S_ell(k, a=2) = 3`, so the average over the units is
  3/2 at every level, where a convex combination would need 1. (An
  earlier version of this bullet said the average IS 1, which
  `main.tex`'s own `rem:no-monotone-certificate` refutes and has since
  2026-08-10; corrected in round 5. See notes/H-022.md, closed.) Admissibility of `j` for fixed `a` depends only on
  `a mod 3` (the condition `2^(j+1)*a = 1 mod 3` fixes the parity of
  `j`, nothing more), so `S_ell(k,a)` takes only two values per level,
  one per residue class of `a mod 3`. Exact counterexample at
  `ell = 3, k = 3`: `S_3(3,a) = 9/7` for `a = 1 mod 3`, `12/7` for
  `a = 2 mod 3` (corrected 2026-08-10: an earlier version listed six
  values, `9/7, 12/7, 6/7, 6/7, 3/7, 6/7` for `a=1,2,4,5,7,8`, i.e.
  four DISTINCT values spread across only two residue classes mod 3,
  which the formula does not permit; the first two are right and the
  rest were wrong. See H-013.md's own 2026-08-10 correction section). H-166's mechanism worked because
  its denominator transforms under the *same* operator as its
  numerator, pointwise; here it does not.
- The diagonal `k = ell + d` is not preserved by the one-step recursion
  (`k - j = (ell-1) + (d+1-j)` mixes many offsets), so any bound this
  route gives for `min_a R_ell(ell+d,a)` degenerates to zero.
- `min_a R_ell` is not even monotone at small levels: exact fractions
  `m_ell = 9/28, 648/1459, 2430/5057, 13851/30695` for `ell = 4..7`
  (these cross-check against committed data: they equal 0.3214,
  0.4441, 0.4805, 0.4512, which is exactly the `d = +12` column of
  `ratio_ell18_output.txt` at those levels. The `S_3(3,a) = 9/7, 12/7`
  values are now produced by `weight_asymmetry.py` in this directory,
  written 2026-08-16 during round 4 after round 3's "documented as
  hand-derived" was correctly judged not to satisfy Rule 9a. It works
  in exact rationals and reproduces 9/7 and 12/7 with Haar mean 1 at
  `ell = 3`, and 6/5 and 9/5 with mean 1 at `ell = 4`. Output in
  `weight_asymmetry_ell3.log`. `notes/H-013.md`'s earlier citation of
  an `exact_weight_asymmetry.py` at a nonexistent path is superseded)
  (matching this project's own floating-point values independently),
  rising then falling.
- What does survive is a genuine convex combination of full probability
  *vectors* across residues, which yields convex-functional bounds
  (entropy, norms), not a positive coordinate minimum.

No route from the composition-counting recursion alone, without new
input, currently closes H-013 either way. See `notes/H-013.md`'s
2026-08-10 section for the full writeup and closure as
`closed-inconclusive`.

## Disk-backed mode for runs that don't fit in RAM

`central_ratio.py` now accepts `--scratch-dir PATH`. When set, the
`(k_max+1) x 3^ell` tables are allocated as disk-backed `np.memmap`
files under `PATH` instead of in RAM; only the current and immediately
previous level's files exist on disk at any time, each deleted as soon
as the next level's file is built and flushed. Without the flag,
behaviour is unchanged from before (arrays in RAM).

A preflight check estimates the peak of the last level transition
(the final level's own table PLUS the still-live previous level's
table, which coexist until `extend()` returns and the old file is
deleted, not the final level's table alone: an undercount of about 33%
found and fixed 2026-08-15 during critique, since the original formula
fed both the RAM-mode `SystemExit` guard below and the disk-backed
free-space check) and, in RAM mode, refuses to start (`SystemExit`,
before allocating anything) if that peak exceeds 90% of currently
available RAM (`--ram-safety-margin` to change the threshold), rather
than risking an `MemoryError` mid-run or locking up the machine.
Disk-backed mode now also refuses to start if the same peak exceeds
free disk space at `--scratch-dir`, rather than only printing an
estimate. Disk-backed mode is I/O-bound rather than memory-bound. What the
`ell = 18` run below actually measured, all within that single
disk-backed run (`ratio_ell18_output.txt`): `ell = 17` took 485.4 s
and `ell = 18` took 3961.0 s, a jump of 8.2x for 3x the data. The
cost is therefore not a flat disk-mode penalty: at `ell = 17` the
table still fits in page cache and disk mode costs little, while at
`ell = 18` it does not and the I/O bottleneck dominates. (An earlier
version of this paragraph compared 3961 s against "536 s for ell=17 in
RAM mode"; that 536 s figure came from a log file that was never
committed to this repository and cannot be reproduced from anything
here, so the like-for-like within-run comparison above replaces it.
Corrected 2026-08-16.)

What disk mode buys is not needing that peak in RAM, about 119 GiB
(89.5 GiB for the `ell=18` table alone, plus the 29.8 GiB `ell=17`
table still open during the transition). **Derived from the model, not
measured**: the scratch directory was deleted after the run and its
peak size was never recorded, so nothing here confirms it
independently.

An earlier version of this sentence claimed the figure "matches the
scratch directory's actual peak disk usage, observed live during the
run". That claim was removed from the evidence block in round 3 and
survived HERE through rounds 3, 4 and 5, being flagged each time. It
survived because the fixes were written against a differently-wrapped
copy of the sentence and silently matched nothing. Fourth occurrence;
recorded because the failure mode is the fix method, not the text.

Verified correctness-preserving on 2026-08-15: `--max-ell 8` in RAM
mode and in `--scratch-dir` mode produce byte-identical output (every
ratio column and every row-sum assertion against
`bounded_composition_counts`), aside from the expected preflight-mode
announcement line. Only the final level's `.dat` file was left behind
by the disk-backed run, as designed.

## Memory and runtime

All figures in this section were recomputed 2026-08-16 against the
corrected peak model (a transition holds the new level's array plus
the still-live previous one, 4/3 of the new array, not two equal
arrays; see the disk-backed section above) and against this
repository's own committed output files, replacing earlier estimates
that used the superseded two-equal-arrays model and pre-migration
timings.

The support run (`central_zeros.py`) holds, during a level transition,
the new level's `3^ell` array of 32-bit words plus the previous
level's `3^(ell-1)` one (so 4/3 of the new array, not two equal
arrays, corrected 2026-08-16), plus a few `3^ell` Boolean arrays; its own `rss=` column records 3.13 GiB at `ell = 18`, 9.09 GiB
at `ell = 19` and 27.11 GiB at `ell = 20`, with those levels taking
99.2 s, 296.2 s and 1010.1 s (`zeros_extended_output.txt`). The
exact-count run (`central_ratio.py`) holds `(k_max+1) x 3^ell` int64
arrays, so the `--offsets 0 5 12 --sqrt-multiples 1` command above, whose
ceiling is `k_max = 28`, peaks at about 12.4 GiB (9.3 GiB final table
plus its predecessor). The other `ell = 16` command listed, the
1458-integer fixed-set one, has `k_max = 21` and a 9.41 GiB peak by the
same model; note that 9.4 GiB is also what `ratio_fixedsets_ell16.log`
prints, but that log is from the 54/486/4374 runs, NOT from this
command, whose own output is carried by `ratio_ell18_output.txt`
(preflight 89.5 GiB, since that run went to `ell = 18`). Round 3's
disambiguation attached the figure to the wrong log and invented a
"third" command; corrected round 4.
`central_quantiles.py` holds the same shape of array, so its memory
scales the same way with `--max-ell` and the largest `--offsets` value.
The `ell = 17` extension of `central_ratio.py` peaks at about 38.5 GiB
(`k_max = 29`), derived from the model rather than from a recorded
measurement (see the note in the ell=17 section above); `ell = 18` (`k_max = 30`, one higher, since `k_max` grows with
`--max-ell`) peaks at about 119 GiB, run disk-backed via
`--scratch-dir` (see below). All three scripts are single-threaded numpy. Lower `--offsets` to
lower `k_max` and the memory falls proportionally.

## `ell = 18` result (`central_ratio.py`, disk-backed)

Run on 2026-08-15 with `--scratch-dir`, extending the same fixed-integer
test one more level. Both tracked columns continue the established
trend with no bend: `d = +12` falls from 0.2378 (`ell = 17`) to 0.2222
(`ell = 18`), a ratio of 0.9344, matching the ~0.94/level decay rate
already documented above (single-step ratios of the `d = +12`
column over `ell = 10..18` vary from 0.907 to 0.993, so this is
unremarkable noise, not a new signal; over the wider `ell = 6..16`
window cited earlier the top of that range is 0.994); the minimum over the 1458 fixed integers at `d = +5`
reads 0.0799, inside that same column's noisy flat band, 0.0736-0.1427
over `ell = 9..18`. (An earlier version of this sentence gave that band
as "0.02-0.14"; the 0.02 end belongs to no column here. The other
`d = +5` column in the same table, the minimum over ALL units, spans
0.0197-0.0536 over `ell = 10..18` and never approaches 0.14.
Corrected 2026-08-16.)
Full row in the evidence block below. Level 18's whole iteration took
3961.0 s (~66 min); that is what the script's per-level timer covers,
which is the `extend()` call plus the level's
`bounded_composition_counts`, its row-sum assertion and its row reads,
not `extend()` alone; levels 1-17 together
took 656.0 s, summing that file's own per-level timings (corrected
2026-08-16 from a stated "about 685 s", which matched no figure in
the file).

## Evidence (Rule 9a): `central_ratio.py --max-ell 18`

```
Command:      python3 central_ratio.py --max-ell 18 --offsets 5 12
              --sqrt-multiples 1 --fixed-offset 5 --random-fixed 1458
              --random-bound 2187 --scratch-dir <scratch>
Commit:       cfc63a1d21ebbfaca9bd93c6787d0715984d2975
Date:         2026-08-15
Environment:  Linux 7.0.0-28-generic x86_64, Python 3.12.3, 62 GiB RAM,
              disk-backed mode (estimated peak ~119 GiB on disk, counting
              the two memmapped tables only: the ell=18 table, 89.5 GiB,
              plus the still-open ell=17 table, 29.8 GiB, during the
              transition. Derived, not measured: the scratch directory was
              deleted after the run and its peak size was not recorded)
Exit:         0
Output:       experiments/E-010-wirsching-conj2-central-zeros/ratio_ell18_output.txt
Checked:      producer (main session), against the raw output file and
              against the established d=+12 decay-rate and fixed-integer
              flat-band pattern from ell=4..17; critiqued twice
              (2026-08-15, 2026-08-16), see notes/H-015.md
```

**The recorded output file predates the peak-memory fix, by design.**
`ratio_ell18_output.txt`'s own first line reads "[preflight] final
level needs ~89.5 GiB on disk", the pre-fix undercount, because the
run executed against commit `cfc63a1` (the Commit field above), before
`a88d674` corrected the formula. Re-running the recorded Command
against HEAD prints "last transition peaks at ~119.3 GiB" on that line
instead. Every other line of the file, all the actual ratio data, is
unaffected: the fix touched only the preflight estimate and the guard
thresholds, never the computation. Noted here rather than regenerating
the file, since a 66-minute rerun would reproduce the same ratio
columns below the first line (the per-level elapsed-time fields vary
between runs, so the file would not be byte-identical).

This closes the `ell = 17` -> `ell = 18` half of H-015's backfill (the
`central_ratio.py` heavy rerun). Formal hypothesis closure still awaits
three consecutive clean critique rounds (Rule 8f); see
`notes/H-015.md` for the running record.


Migration note: references `H-166`, out of this repo's scope, left unrenumbered.

## Evidence (Rule 9a): `weight_asymmetry.py`

```
Command:      python3 weight_asymmetry.py --ell 3 --k 3
              python3 weight_asymmetry.py --ell 4 --k 4
Commit:       1e30233
Date:         2026-08-16
Environment:  Linux 7.0.0-28-generic x86_64, Python 3.12.3; exact
              rationals via fractions.Fraction, no floating point
Exit:         0 (both)
Output:       experiments/E-010-wirsching-conj2-central-zeros/weight_asymmetry_ell3.log
              experiments/E-010-wirsching-conj2-central-zeros/weight_asymmetry_ell4.log
Checked:      reproduces the 9/7 and 12/7 values this document and
              notes/H-013.md have quoted since 2026-08-10, and gives
              6/5, 9/5 at ell=4. Critiqued in H-015 round 5, which found
              the first version's Haar-mean claim wrong and a gbar edge
              case at ell=0 that broke the identity at ell=1. The script
              now ASSERTS main.tex's identity S(a=1)+S(a=2)=3 and the
              unit mean 3/2, and has been run at ell=1..7.
```

**What this script establishes, after two corrections.** `notes/H-013.md` writes the weight with a leading
factor of 3 and asserts the weight sum averages to 1. Computed exactly:
with the factor the values are 27/7 and 36/7, without it 9/7 and 12/7,
which is what the same file quotes two lines below its own formula. So
formula and quoted values already disagree.

The factor turns out to be **correct** under the reading that `bar_g`
means the raw count `B` there, since `B_ell = 2*3^(ell-1)*gbar_ell`
makes `3*B_(ell-1)/B_ell` identically `gbar_(ell-1)/gbar_ell`. The
notation is overloaded, not wrong.

The error was the separate claim that the average is 1. It is 3/2 over
the units, by the identity `S(a=1) + S(a=2) = 3` that `main.tex`'s
`rem:no-monotone-certificate` proves from the bounded-urn recursion.
This document and the script's first version both reported "mean 1",
which was the average over ALL residues: the inadmissible class
contributes an artificial zero and the units are exactly 2/3 of the
total, so that average is always 2/3 of the unit average. Obtaining 1
that way follows from the measure chosen, not from an identity.

Tracked and resolved as **H-022** (opened and closed 2026-08-16). H-013's
closure conclusion is unaffected: what that argument needs is that the
weight sum is not constant in `a`, which the identity establishes.
