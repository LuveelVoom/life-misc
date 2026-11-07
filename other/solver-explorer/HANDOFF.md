# Solver run explorer — chat handoff

Upload this whole folder (or the zip) to a new chat and paste the "Starter prompt" below.

## Files
| File | What it is |
|---|---|
| `solver_log.jsonl` | The original solver log you pasted (492 JSONL events, 2026-09-19 18:16 → 2026-09-25 01:09 UTC). |
| `solver_run_explorer.html` | The finished tool, latest version (includes the interruption filter). Self-contained: open in any browser. Needs internet only to load d3 v7.8.5 from cdnjs. |
| `HANDOFF.md` | This file. |

Note: the HTML embeds the log in a compact hand-transcribed form (`<script type="text/plain" id="data">`, lines like `R|19T18:16:55|root|54|memory_cap`). For exactness, use the page's "Load another log (.jsonl)" button with `solver_log.jsonl`. The compact format drops the 480-byte final-stub checkpoint deletions.

## What was asked and built
1. **Request:** "a thorough graphing tool element of various factors (partitions, time between partitions, partition graph, depths, depths of each branch, etc.)" from the log.
2. **Follow-up:** "add an option to hide things that went overly long due to interrupts" — a ~68 h interrupt was squashing the time charts.

### Views in the tool (all linked; click a branch to highlight it everywhere)
- Summary stats strip
- Search tree by height (branch start→end height; partitions as joins; queued branches as open circles)
- Activity timeline (Extend lane, Partition lane, manager-pause triangles)
- Height over time (run start→end segments; straight line is an assumption)
- Depth of each branch (bars, in run order)
- Partition build time (split request → materialized)
- Time between partitions (labelled with split height)
- Run length vs levels gained (scatter)
- Checkpoint sizes over time (split/resume vs partition base)
- Sortable partition and branch tables
- "Load another log (.jsonl)" button

### Interruption filter (second version)
- Checkbox "Hide time inflated by interruptions and pauses" (default on). Hides runs resumed after an `interrupted` result (resume time isn't logged) and partition builds that included a `manager_paused` event; flags/hides partition gaps spanning either.
- "Collapse idle gaps longer than N h" (default 2): time axes become piecewise, with each long empty stretch shrunk to a grey ≈ band (hover shows hours hidden).
- Each chart notes how many items were hidden.
- Tested only on a synthetic copy of the log with a 66 h delay inserted. The user's real 68 h-interrupt log was NOT seen; if it logs resumes differently (e.g. a `branch_started` on resume or a new event kind), the filter may need adjusting. Open item: ask for a few lines around that interrupt.

## How the log was interpreted (parsing rules)
- Event kinds: `solver_result` (reason: `memory_cap`, `exhausted`, `interrupted`, `row_limit`; stage `extend` or `materialize`), `branch_split_requested`, `manager_paused` (stage `partition` or `extend`), `partition_materialized` (always 4 children), `branch_started`, `checkpoint_deleted` (bytes, source).
- A branch hitting `memory_cap` splits into `partition-N/1..4`, children start at split height + 1. Branches run sequentially, depth-first.
- `materialize`/`row_limit` results are child creation, not runs.
- After `interrupted`, the next result for that branch has no logged start; drawn from the interruption time.
- Root's start height isn't logged; drawn from 0 and shown as "before log".
- Checkpoint classes: `partitions/...` = partition base; `checkpoints/...` = split/resume state; 480-byte = final stub.

## Key results from this log (computed, high confidence)
- 126.9 h wall time; ~108 h extending, ~19 h building partitions (raw, before filtering).
- 23 partitions; 63 branches exhausted; 23 split; 0 completions.
- Max height 659 (branch 5/1).
- 7 branches still queued: 1/3, 1/4 (start 55); 12/3, 12/4 (61); 14/3, 14/4 (72); 23/4 (158).
- Typical partition build 17–56 min; outliers 1.9 h, 2.8 h, 5.0 h all included a manager pause.
- ~0.52 TB of checkpoints deleted in total.

## Your stated chat preferences (set these on the new account under Settings → Profile if you want them applied)
Be concise, no repetition; no unfounded "deep" generalizations; don't agree reflexively when corrected — assess both sides; treat each chat as independent (no cross-chat context); state confidence (high / moderate / low) and source types.

## Starter prompt
> I'm continuing work from another chat. Attached: HANDOFF.md (context), solver_log.jsonl (original data), solver_run_explorer.html (the current tool). Read HANDOFF.md first, then publish solver_run_explorer.html as an artifact unchanged so I can keep iterating on it. Next I want to: [your next request].
