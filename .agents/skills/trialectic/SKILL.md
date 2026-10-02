---
name: trialectic
description: Use for explicit GRHayLET trialectic, three-seat independent review, "use tri", "run tri", or "engage tri", or when repository policy requires it. Mentioning, auditing, or editing this skill does not activate it. Three independent seats; delta-and-impact follow-ups; verified /work/ delivery.
---

# GRHayLET Trialectic

Read [the shared protocol](../../review-protocol.md) with this skill; it owns modes,
independence, coverage, budgets, validation, and delivery. Distribute both and the
[GRHayLET reference](../../grhaylet-evidence.md); load only relevant reference sections,
not the other skill.

Use exactly three independent seats, excluding root. Root owns live writes and
delivery; seats never delegate or edit live targets. Default to `review`, including
documents and all KB work. Use `design` or `tri-build` only under protocol criteria.

## Roles

**Physics, numerics, and behavior.** Check relevant GRHD/GRMHD, initial-data, EOS,
leakage, reconstruction, and induction assumptions against the owning thorn's source.
Check units, indices, signs, densitization, face/stagger ownership, bounds, recovery
paths, stability, tolerances, and scientific oracles for the applicable driver/mode.

**API, build, and integration.** Check thorn headers, `interface.ccl`, `param.ccl`,
`configuration.ccl`, recursive `make.code.defn`, scheduling/storage/synchronization,
Carpet/CarpetX differences, GRHayLib call sites, and HydroBase/ADM/MoL/Tmunu boundaries.
Prefer existing GRHayLET machinery and the simplest sufficient change; external
GRHayL behavior needs separate proof.

**Tests, docs, and evidence.** Trace criteria to Cactus test declarations, parfiles,
checked-in observations, tolerances, thorn documentation, and branch-owned KB routes.
Check source registration, stable locators, reverse dependencies, issues, and the
canonical KB linter where applicable. Catch coverage/delivery omissions and claims
that outrun static or runtime evidence; do not audit unrelated repository areas.

For non-code work, use a first-principles analyst, a domain/implementation expert,
and an adversarial evidence editor. Roles add complementary scrutiny, not duplicate
suites.

## Completion

Follow-ups cover deltas, open findings, and affected code/docs, retaining valid
coverage. Root must [deliver and verify](../../review-protocol.md#delivery) ALL
intended authorized `/work/` destinations. A verdict, `DRAFT COMPLETE`, or exhausted
budget is not delivery. Install acceptable work; preserve blocked drafts separately.

If the protocol is missing, preserve drafts/companions in fresh
`/work/review-results/<task>/unapproved/` outside active source/KB/instruction discovery;
verify copies and report the missing file, never invent approval. Explicit outputs/
scoped write limits control; report inaccessible storage and actual retained paths.
