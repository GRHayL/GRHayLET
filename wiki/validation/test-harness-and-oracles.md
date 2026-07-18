# Test Harness and Oracles

> Status: confirmed · Last reconciled: 07-17-2026
> Up: [Validation](index.md)

## Summary

IllinoisGRMHD ships one `test.ccl`, six top-level test parameter files, five
checked-in per-case parameter fixtures, and 40 `.asc` files. They define
configured comparisons and generated reference evidence. No in-scope build
workflow or latest test result exists, so file presence is not proof of a
current pass.

Claim evidence:

- Claim: Checked-in `.asc` files are generated evidence for visible fixture
  context, not proof of a current test pass; no in-scope current result exists.
- Role: generated evidence
- Deciding authority: `IllinoisGRMHD/test/test.ccl` comparison declarations
  and `IllinoisGRMHD/test/**/*.asc` oracle roles and headers
- Corroboration: no checked-in current command/environment/result record is available
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-run; tool_version=not-run; backend=not-run; precision=not-run; GPU=not-run; restart=not-run; distributed=not-run; error_path=not-run; options=static inventory only; date=07-17-2026`

## Detail

### `test.ccl` contract

Thorn-wide tolerances are `ABSTOL 1e-11` and `RELTOL 1e-11`.

| Declaration | State | Process setting |
| --- | --- | --- |
| `TEST TOV` | active | no test-specific `NPROCS` |
| `TEST Balsara1` | active | `NPROCS 1` |
| `TEST Balsara2` | active | `NPROCS 1` |
| `TEST Balsara3` | active | `NPROCS 1` |
| `TEST Balsara4` | commented out | commented `NPROCS 1` |
| `TEST Balsara5` | active | `NPROCS 1` |

Comment above Balsara4 says test is too sensitive while using PPM
reconstruction. This is checked-in rationale, not measured diagnosis from a
new run.

Official Cactus Users' Guide says `test.ccl` can set thorn/test tolerances and
processor count, `TEST` names test examples, and reference ASCII output is
placed under directory named for parameter-file base. It also says absent
`NPROCS` means a test is assumed to run with any processor count. Those are
generic harness semantics only. Repository-specific active rows and files
come from local tree.

### File-role inventory

| Layer | Count | Role |
| --- | ---: | --- |
| `par/*.par` | 6 | Shipped example configurations: Balsara1–5 and `magnetizedTOV`. Not harness declarations. |
| `test/*.par` | 6 | Top-level checked-in test inputs/wrappers: Balsara1–5 and `magnetizedTOV`. |
| `test/<case>/<case>.par` | 5 | Checked-in per-case fixtures for Balsara1, 2, 3, 5, and `magnetizedTOV`; their headers report automatic generation by Cactus, but this KB does not infer which harness action created them. |
| `test/<case>/*.asc` | 40 | Checked-in generated evidence used as reference outputs/oracles for visible fixture context. |
| `test/test.ccl` | 1 | Test registration and comparison options. |

Total visible inventory: 52 files. Balsara1, 2, 3, and 5 each have eight
profile files: `rho`, `press`, `vx`, `vy`, `vz`, and centered `Bx`, `By`, `Bz`,
all `.x.asc`. Magnetized TOV has minimum and maximum scalar files for `rho`
and each centered B component, eight files total.

`.asc` headers identify CarpetIOASCII or CarpetIOScalar, source parameter
filename, variables, and recorded samples. They are generated evidence for
that historical context. Without executing comparison in a recorded build,
they establish neither current pass nor numerical validity. KB work must not
regenerate or overwrite them.

### Visible gaps

- Balsara4 has example and top-level test `.par` files, but its `TEST` block
  is commented and no `test/Balsara4/` oracle directory is visible.
- `test.ccl` declares `TEST TOV`, while visible top-level file and fixture
  directory are named `magnetizedTOV`. This tree does not define Cactus test
  discovery/name mapping, so this is naming uncertainty—not a claimed harness
  failure or contradiction. Resolution needs authorized harness enumeration
  in configured Cactus environment.
- No checked-in command log, environment record, revision-specific result, or
  current build output is in scope. Thus no relation is labeled “Validated by.”

## Sources

- [`IllinoisGRMHD/test/test.ccl`](../../IllinoisGRMHD/test/test.ccl) — exact
  tolerances, active/commented `TEST` blocks, process settings, and Balsara4 comment.
- [`IllinoisGRMHD/test/`](../../IllinoisGRMHD/test/) — complete checked-in
  test input, fixture, and oracle path inventory.
- [`IllinoisGRMHD/test/Balsara1/rho.x.asc`](../../IllinoisGRMHD/test/Balsara1/rho.x.asc) —
  representative CarpetIOASCII generated-evidence header and samples.
- [`IllinoisGRMHD/test/magnetizedTOV/rho.maximum.asc`](../../IllinoisGRMHD/test/magnetizedTOV/rho.maximum.asc) —
  representative CarpetIOScalar generated-evidence header and samples.
- [Cactus Users' Guide, “Adding a Test Suite”](https://einsteintoolkit.org/usersguide/UsersGuide.html#x1-350000C1.8.5) —
  official generic `test.ccl`, tolerance, `NPROCS`, and reference-output semantics.

## See Also

- Parent: [Validation](index.md)
- Example: [Balsara and TOV Cases](balsara-and-tov-cases.md)
- Depends on: [Schedule Lifecycle](../architecture/schedule-lifecycle.md)
- See also: [Reconstruction, Fluxes, and Sources](../evolution/reconstruction-fluxes-and-sources.md)
