// transcensus-int: port of transcensus.c (orig/) to Life-like, isotropic non-totalistic (Hensel) and isotropic MAP rules
// The census itself is unchanged: same soups (same random stream), same grid and edge handling, same state machine, same canonical forms
// and the same LifeHistory output. With the default rule B3/S23 and default options the output matches the original exactly.
// Additions: a rule option, command-line versions of the constants, deterministic multithreading and a plain RLE list output

#include <stdlib.h>
#include <inttypes.h>
#include <memory.h>
#include <time.h>
#include <string.h>
#include <stdio.h>

#ifdef _OPENMP
	#include <omp.h>
#endif

#include "lib/lib.c"
#include "lib/rect.c"
#include "lib/celllist.c"
#include "lib/objects.c"
#include "lib/randomarray.c"
#include "lib/rule.c"
#include "lib/golgrid.c"
#include "lib/gridmisc.c"
#include "lib/hashtable.c"
#include "lib/golutils.c"
#include "lib/goltext.c"

#define GRID_SIZE 192
#define MAX_MIN_PATTERNS 16
#define MAX_STATS_POP 96
#define MAX_STATE_A_POP_THRESHOLD 12
#define CHUNK_SOUPS 4096

// Defaults are the constants of the original
static s64 used_soup_cnt = 10000000;
static s64 max_key_patterns = 2000000;
static s64 report_threshold = 200;
static s32 state_a_pop_threshold = 12;
static s32 state_b_pop_threshold = 20;
static s32 state_c_unique_pop_threshold = 4;
static s32 max_gens = 256;

typedef struct
{
	s64 cnt;
	u64 hash;
	ObjCellList_Cell min_pattern_cell [MAX_STATE_A_POP_THRESHOLD];
	ObjCellList min_pattern;
} MinPattern;

typedef struct
{
	s64 cnt;
	s32 min_pattern_cnt;
	MinPattern min_pattern [MAX_MIN_PATTERNS];
} KeyPattern;

// What a worker thread hands to the in-order merge for one settled soup
typedef struct
{
	s32 soup_offset;
	u64 key_hash;
	u64 min_hash;
	Rect min_rect;
	s32 min_cell_cnt;
	ObjCellList_Cell min_cell [MAX_STATE_A_POP_THRESHOLD];
} SettledRecord;

typedef struct
{
	GoLGrid *p_m1;
	GoLGrid *p_m0;
	GoLGrid *soup;
	GoLGrid *key_pattern;
	GoLGrid *min_pattern;
	GoLGrid *canonical_key_pattern;
	GoLGrid *canonical_min_pattern;
	GoLGrid *canonical_t1;
	GoLGrid *canonical_t2;
} ThreadContext;

typedef enum {STATE_PRE_A, STATE_A, STATE_B, STATE_C} ObjectIdentificationState;

static const s64 min_value_multiplier [MAX_STATE_A_POP_THRESHOLD + 1] = {128, 128, 128, 128, 128, 96, 48, 24, 12, 7, 5, 4, 3};


static __not_inline GoLGrid *alloc_grid (void)
{
	return GoLUtils_alloc_std_grid (-(GRID_SIZE / 2), -(GRID_SIZE / 2), GRID_SIZE, GRID_SIZE);
}

static __not_inline int run_pattern (ThreadContext *ctx, const GoLGrid *start_pattern, GoLGrid *key_pattern, GoLGrid *min_pattern)
{
	GoLGrid_copy_noinline (start_pattern, ctx->p_m0);

	s32 pop_stats [MAX_STATS_POP];
	s32 unique_pops = 0;

	s32 pop_ix;
	for (pop_ix = 0; pop_ix < MAX_STATS_POP; pop_ix++)
		pop_stats [pop_ix] = 0;

	s32 gen = 0;
	s32 min_pop = -1;
	ObjectIdentificationState state = STATE_PRE_A;

	while (TRUE)
	{
		s32 pop = GoLGrid_get_population (ctx->p_m0);

		if (state == STATE_PRE_A)
			if (pop <= state_a_pop_threshold)
				state = STATE_A;

		if (state == STATE_A)
		{
			if (min_pop == -1 || pop < min_pop)
			{
				min_pop = pop;
				GoLGrid_copy (ctx->p_m0, min_pattern);
			}

			if (pop >= state_b_pop_threshold)
			{
				state = STATE_B;

				if (unique_pops > 0)
					for (pop_ix = 0; pop_ix < MAX_STATS_POP; pop_ix++)
						pop_stats [pop_ix] = 0;

				unique_pops = 0;
			}
		}
		else if (state == STATE_B)
		{
			if (pop <= state_a_pop_threshold)
				state = STATE_PRE_A;
			else
				if (pop < MAX_STATS_POP)
				{
					if (pop_stats [pop] == 0)
						unique_pops++;

					pop_stats [pop]++;

					if (unique_pops >= state_c_unique_pop_threshold)
					{
						GoLGrid_copy (ctx->p_m0, key_pattern);
						return TRUE;
					}
				}
		}

		if (gen >= max_gens)
			return FALSE;

		GoLGrid *temp = ctx->p_m1;
		ctx->p_m1 = ctx->p_m0;
		ctx->p_m0 = temp;

		GoLGrid_evolve (ctx->p_m1, ctx->p_m0);
		gen++;
	}
}

static __not_inline int remove_oldest_single (KeyPattern *key_pattern)
{
	s32 min_ix;
	s32 src_ix;
	for (min_ix = 0; min_ix < key_pattern->min_pattern_cnt; min_ix++)
		if (key_pattern->min_pattern [min_ix].cnt == 1)
		{
			for (src_ix = min_ix + 1; src_ix < key_pattern->min_pattern_cnt; src_ix++)
			{
				key_pattern->min_pattern [src_ix - 1].cnt = key_pattern->min_pattern [src_ix].cnt;
				key_pattern->min_pattern [src_ix - 1].hash = key_pattern->min_pattern [src_ix].hash;
				ObjCellList_copy (&key_pattern->min_pattern [src_ix].min_pattern, &key_pattern->min_pattern [src_ix - 1].min_pattern);
			}

			key_pattern->min_pattern_cnt--;
			return TRUE;
		}

	return FALSE;
}

// Same as the original store_min_pattern, except that the cell list was already made from the canonical min pattern by the worker thread
static __not_inline void store_min_pattern (const SettledRecord *rec, KeyPattern *key_pattern)
{
	s32 min_ix;
	for (min_ix = 0; min_ix < key_pattern->min_pattern_cnt; min_ix++)
		if (key_pattern->min_pattern [min_ix].hash == rec->min_hash)
		{
			key_pattern->min_pattern [min_ix].cnt++;
			return;
		}

	if (key_pattern->min_pattern_cnt >= MAX_MIN_PATTERNS)
		if (!remove_oldest_single (key_pattern))
			return;

	MinPattern *mp = &key_pattern->min_pattern [key_pattern->min_pattern_cnt];
	mp->cnt = 1;
	mp->hash = rec->min_hash;

	ObjCellList_make_empty (&mp->min_pattern, mp->min_pattern_cell, MAX_STATE_A_POP_THRESHOLD);
	s32 cell_ix;
	for (cell_ix = 0; cell_ix < rec->min_cell_cnt; cell_ix++)
		mp->min_pattern_cell [cell_ix] = rec->min_cell [cell_ix];
	mp->min_pattern.cell_cnt = rec->min_cell_cnt;
	Rect_copy (&rec->min_rect, &mp->min_pattern.obj_rect);

	key_pattern->min_pattern_cnt++;
}

static __not_inline s64 get_min_pattern_value (const MinPattern *mp)
{
	s32 pop = mp->min_pattern.cell_cnt;
	return mp->cnt * min_value_multiplier [pop];
}

static __not_inline int compare_key_patterns (const void *key_pattern_1, const void *key_pattern_2)
{
	if ((*(const KeyPattern **) key_pattern_1)->cnt < (*(const KeyPattern **) key_pattern_2)->cnt)
		return 1;
	else if ((*(const KeyPattern **) key_pattern_1)->cnt > (*(const KeyPattern **) key_pattern_2)->cnt)
		return -1;
	else
		return 0;
}

static __not_inline void create_sorted_key_patterns (s32 key_pattern_cnt, KeyPattern *key_pattern_store, KeyPattern **sorted_key_pattern)
{
	s32 key_ix;
	for (key_ix = 0; key_ix < key_pattern_cnt; key_ix++)
		sorted_key_pattern [key_ix] = &key_pattern_store [key_ix];

	fprintf (stderr, "Sorting...\n");
	qsort (sorted_key_pattern, key_pattern_cnt, sizeof (sorted_key_pattern [0]), &compare_key_patterns);
	fprintf (stderr, "Sorted\n");
}

static __not_inline int ThreadContext_create (ThreadContext *ctx)
{
	ctx->p_m1 = alloc_grid ();
	ctx->p_m0 = alloc_grid ();
	ctx->soup = alloc_grid ();
	ctx->key_pattern = alloc_grid ();
	ctx->min_pattern = alloc_grid ();
	ctx->canonical_key_pattern = alloc_grid ();
	ctx->canonical_min_pattern = alloc_grid ();
	ctx->canonical_t1 = alloc_grid ();
	ctx->canonical_t2 = alloc_grid ();

	return ctx->p_m1 && ctx->p_m0 && ctx->soup && ctx->key_pattern && ctx->min_pattern && ctx->canonical_key_pattern && ctx->canonical_min_pattern &&
			ctx->canonical_t1 && ctx->canonical_t2;
}

// Runs soup_cnt soups and writes one record per settled soup, in soup order. Returns the number of records
static __not_inline s32 run_chunk (ThreadContext *ctx, const RandomDataArray *rda, const u64 *soup_bits, s32 soup_cnt, SettledRecord *rec)
{
	s32 rec_cnt = 0;
	s32 soup_offset;
	for (soup_offset = 0; soup_offset < soup_cnt; soup_offset++)
	{
		GoLGrid_clear_noinline (ctx->soup);
		GoLGrid_or_8_by_8_block_noinline (ctx->soup, 0, 0, soup_bits [soup_offset]);

		if (!run_pattern (ctx, ctx->soup, ctx->key_pattern, ctx->min_pattern))
			continue;

		SettledRecord *r = &rec [rec_cnt++];
		r->soup_offset = soup_offset;
		r->key_hash = GoLUtils_make_canonical (ctx->key_pattern, ctx->canonical_key_pattern, rda, ctx->canonical_t1, ctx->canonical_t2);
		r->min_hash = GoLUtils_make_canonical (ctx->min_pattern, ctx->canonical_min_pattern, rda, ctx->canonical_t1, ctx->canonical_t2);

		ObjCellList ocl;
		ObjCellList_make_empty (&ocl, r->min_cell, MAX_STATE_A_POP_THRESHOLD);
		GoLGrid_to_obj_cell_list (ctx->canonical_min_pattern, &ocl);
		r->min_cell_cnt = ocl.cell_cnt;
		Rect_copy (&ocl.obj_rect, &r->min_rect);
	}

	return rec_cnt;
}

// One-line RLE of a cell list in the census rule, for the --list output
static __not_inline void print_obj_cell_list_rle (FILE *stream, const ObjCellList *ocl)
{
	char grid [GRID_SIZE] [GRID_SIZE + 1];
	s32 width = 0;
	s32 height = 0;
	s32 cell_ix;
	for (cell_ix = 0; cell_ix < ocl->cell_cnt; cell_ix++)
	{
		width = higher_of_s32 (width, ocl->cell [cell_ix].x + 1);
		height = higher_of_s32 (height, ocl->cell [cell_ix].y + 1);
	}

	s32 y;
	s32 x;
	for (y = 0; y < height; y++)
		for (x = 0; x < width; x++)
			grid [y] [x] = 'b';
	for (cell_ix = 0; cell_ix < ocl->cell_cnt; cell_ix++)
		grid [ocl->cell [cell_ix].y] [ocl->cell [cell_ix].x] = 'o';

	for (y = 0; y < height; y++)
	{
		s32 end = width;
		while (end > 0 && grid [y] [end - 1] == 'b')
			end--;

		x = 0;
		while (x < end)
		{
			s32 run = 1;
			while (x + run < end && grid [y] [x + run] == grid [y] [x])
				run++;
			if (run > 1)
				fprintf (stream, "%d", run);
			fprintf (stream, "%c", grid [y] [x]);
			x += run;
		}

		if (y < height - 1)
			fprintf (stream, "$");
	}
	fprintf (stream, "!");
}

static __not_inline void print_usage (void)
{
	fprintf (stderr,
		"Usage: transcensus-int [options]\n"
		"  -r RULE   rule: Life-like, Hensel (INT) or isotropic MAP, Golly syntax (default B3/S23). B0 rules are rejected\n"
		"  -n N      settled soups to use (default 10000000)\n"
		"  -t N      report key patterns seen at least N times (default 200)\n"
		"  -g N      give up on a soup after N generations (default 256)\n"
		"  -a N      population at or below which the min pattern is tracked, 1..12 (default 12)\n"
		"  -b N      population at which the pattern counts as active again (default 20)\n"
		"  -c N      distinct populations below 96 needed before the key pattern is taken (default 4)\n"
		"  -k N      key pattern store size; the run stops early and reports if it fills (default 2000000)\n"
		"  -s SEED   seed the soup generator (default: the original fixed state)\n"
		"  -j N      threads (default 1). Output does not depend on the thread count\n"
		"  -l FILE   also write a plain list: count, then the chosen min pattern as RLE in the census rule\n"
		"  --rule-info  print the parsed rule (MAP string, BDD size) and exit\n"
		"  --evolve G   test mode: read an RLE from stdin, place its top-left cell at (0, 0) of the census grid,\n"
		"               run G generations with the census engine and print the result as cells \"x y\" (any rule without B0)\n"
		"  --emit-c     print the rule as C for a build fixed to that rule (see the Makefile) and exit\n"
		"Output (stdout) is the original's: one line per reported key pattern, then a LifeHistory RLE of the min patterns\n"
		"with their counts written above them in state 4. The LifeHistory rule there is only for display.\n");
}

static __not_inline int parse_s64_arg (const char *s, s64 min, s64 max, s64 *out)
{
	u64 n;
	if (!s || !str_to_u64 (s, &n) || n < (u64) min || n > (u64) max)
		return FALSE;
	*out = (s64) n;
	return TRUE;
}

static __not_inline u64 splitmix64 (u64 *x)
{
	u64 z = (*x += 0x9e3779b97f4a7c15u);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9u;
	z = (z ^ (z >> 27)) * 0x94d049bb133111ebu;
	return z ^ (z >> 31);
}

static __not_inline int main_do (s32 thread_cnt, const char *list_file_name)
{
	s32 chunks_per_batch = (thread_cnt == 1 ? 1 : 8 * thread_cnt);

	ThreadContext *ctx = malloc (thread_cnt * sizeof (ThreadContext));
	SettledRecord **chunk_rec = malloc (chunks_per_batch * sizeof (SettledRecord *));
	s32 *chunk_rec_cnt = malloc (chunks_per_batch * sizeof (s32));
	u64 *soup_bits = malloc ((u64) chunks_per_batch * CHUNK_SOUPS * sizeof (u64));
	if (!ctx || !chunk_rec || !chunk_rec_cnt || !soup_bits)
	{
		fprintf (stderr, "Out of memory\n");
		return FALSE;
	}

	s32 thread_ix;
	for (thread_ix = 0; thread_ix < thread_cnt; thread_ix++)
		if (!ThreadContext_create (&ctx [thread_ix]))
			return FALSE;

	s32 chunk_ix;
	for (chunk_ix = 0; chunk_ix < chunks_per_batch; chunk_ix++)
		if (!(chunk_rec [chunk_ix] = malloc (CHUNK_SOUPS * sizeof (SettledRecord))))
		{
			fprintf (stderr, "Out of memory\n");
			return FALSE;
		}

	s32 key_pattern_cnt = 0;
	KeyPattern *key_pattern_store = malloc (max_key_patterns * sizeof (KeyPattern));
	KeyPattern **sorted_key_pattern = malloc (max_key_patterns * sizeof (KeyPattern *));
	if (!key_pattern_store || !sorted_key_pattern)
	{
		fprintf (stderr, "Out of memory for %" PRIi64 " key patterns\n", max_key_patterns);
		return FALSE;
	}

	RandomDataArray rda;
	RandomDataArray_create (&rda, GRID_SIZE * GRID_SIZE / 64);

	HashTable_u64 seen_key_patterns;
	HashTable_u64_create (&seen_key_patterns, 64, 0.75, 0.95);

	s64 settled_cnt = 0;
	s64 soup_ix = 0;
	int done = FALSE;
	int overflow = FALSE;

	while (!done)
	{
		// Soups are drawn from the generator in order, exactly as the original draws one per tested soup
		s32 batch_soups = chunks_per_batch * CHUNK_SOUPS;
		s32 bit_ix;
		for (bit_ix = 0; bit_ix < batch_soups; bit_ix++)
			soup_bits [bit_ix] = random_u64 ();

		#pragma omp parallel for schedule(dynamic, 1) num_threads(thread_cnt)
		for (chunk_ix = 0; chunk_ix < chunks_per_batch; chunk_ix++)
		{
			s32 my_thread = 0;
			#ifdef _OPENMP
				my_thread = omp_get_thread_num ();
			#endif
			chunk_rec_cnt [chunk_ix] = run_chunk (&ctx [my_thread], &rda, soup_bits + (u64) chunk_ix * CHUNK_SOUPS, CHUNK_SOUPS, chunk_rec [chunk_ix]);
		}

		// Merge in soup order, with the original's loop structure, so the result does not depend on the thread count
		for (chunk_ix = 0; chunk_ix < chunks_per_batch && !done; chunk_ix++)
		{
			s32 rec_ix = 0;
			s32 soup_offset;
			for (soup_offset = 0; soup_offset < CHUNK_SOUPS; soup_offset++)
			{
				if (settled_cnt >= used_soup_cnt)
				{
					done = TRUE;
					break;
				}

				if (soup_ix % 100000 == 0)
					fprintf (stderr, "Testing soup %" PRIi64 ", %" PRIi64 " used, %d unique key patterns\n", soup_ix, settled_cnt, key_pattern_cnt);

				soup_ix++;

				if (rec_ix >= chunk_rec_cnt [chunk_ix] || chunk_rec [chunk_ix] [rec_ix].soup_offset != soup_offset)
					continue;

				const SettledRecord *rec = &chunk_rec [chunk_ix] [rec_ix++];

				u64 key_pattern_ix;
				if (HashTable_u64_get_data (&seen_key_patterns, rec->key_hash, &key_pattern_ix))
					key_pattern_store [key_pattern_ix].cnt++;
				else
				{
					if (key_pattern_cnt >= max_key_patterns)
					{
						fprintf (stderr, "Key pattern store overflow: stopping at soup %" PRIi64 " and reporting what was found so far (use -k to raise the limit)\n", soup_ix);
						soup_ix--;
						overflow = TRUE;
						done = TRUE;
						break;
					}

					key_pattern_ix = key_pattern_cnt;
					key_pattern_cnt++;

					HashTable_u64_store (&seen_key_patterns, rec->key_hash, key_pattern_ix, TRUE, NULL);

					key_pattern_store [key_pattern_ix].cnt = 1;
					key_pattern_store [key_pattern_ix].min_pattern_cnt = 0;
				}

				settled_cnt++;
				store_min_pattern (rec, &key_pattern_store [key_pattern_ix]);
			}
		}
	}

	create_sorted_key_patterns (key_pattern_cnt, key_pattern_store, sorted_key_pattern);

	FILE *list_file = NULL;
	if (list_file_name)
	{
		list_file = fopen (list_file_name, "w");
		if (!list_file)
			fprintf (stderr, "Cannot write %s\n", list_file_name);
	}

	// The original's output grid is 16384 high, which holds 24 * 170 reported patterns. Grow it only when more are reported
	s32 report_cnt = 0;
	s32 store_ix;
	for (store_ix = 0; store_ix < key_pattern_cnt; store_ix++)
		if (sorted_key_pattern [store_ix]->cnt >= report_threshold)
			report_cnt++;
	s32 out_height = higher_of_s32 (16384, align_up_s32 (64 + 96 * ((report_cnt + 23) / 24), GOLGRID_HEIGHT_GRANULARITY));

	GoLGrid *out_grid_p = GoLUtils_alloc_std_grid (0, 0, 2048, out_height);
	GoLGrid *out_grid_t = GoLUtils_alloc_std_grid (0, 0, 2048, out_height);

	s32 x_pos = 64;
	s32 y_pos = 64;

	s32 min_ix;
	for (store_ix = 0; store_ix < key_pattern_cnt; store_ix++)
		if (sorted_key_pattern [store_ix]->cnt >= report_threshold)
		{
			printf ("Key pattern %d, count = %" PRIi64 ":\n", store_ix, sorted_key_pattern [store_ix]->cnt);

			s32 max_value_ix = -1;
			s64 max_value = -1;
			for (min_ix = 0; min_ix < sorted_key_pattern [store_ix]->min_pattern_cnt; min_ix++)
			{
				s64 min_pattern_value = get_min_pattern_value (&sorted_key_pattern [store_ix]->min_pattern [min_ix]);
				if (min_pattern_value > max_value)
				{
					max_value_ix = min_ix;
					max_value = min_pattern_value;
				}
			}

			GoLGrid_or_obj_cell_list (out_grid_p, &sorted_key_pattern [store_ix]->min_pattern [max_value_ix].min_pattern, x_pos, y_pos);

			if (list_file)
			{
				fprintf (list_file, "%" PRIi64 "\t", sorted_key_pattern [store_ix]->cnt);
				print_obj_cell_list_rle (list_file, &sorted_key_pattern [store_ix]->min_pattern [max_value_ix].min_pattern);
				fprintf (list_file, "\n");
			}

			char buf [24];
			sprintf (buf, "%" PRIi64, sorted_key_pattern [store_ix]->cnt);
			GoLText_or_text (out_grid_t, buf, x_pos - (3 * strlen (buf)), y_pos - 32, FALSE);

			x_pos += 80;
			if (x_pos >= (2048 - 64))
			{
				x_pos = 64;
				y_pos += 96;
			}
		}

	GoLGrid_print_life_history_full (stdout, NULL, out_grid_p, out_grid_t, NULL, NULL);

	if (list_file)
		fclose (list_file);

	fprintf (stderr, "Used store positions = %d\n", key_pattern_cnt);
	fprintf (stderr, "Settled: %" PRIi64 " of %" PRIi64 "\n", settled_cnt, soup_ix);
	if (overflow)
		fprintf (stderr, "Note: stopped early because the key pattern store was full\n");

	return TRUE;
}

// Test mode for the engine: RLE on stdin (header and # lines skipped, b/. dead, any other letter live), G generations on the
// census grid (192 by 192, dead outside, top-left corner at (-96, -96)), then one "x y" line per live cell
static __not_inline int evolve_test (s64 gens)
{
	GoLGrid *a = alloc_grid ();
	GoLGrid *b = alloc_grid ();
	if (!a || !b)
		return FALSE;
	GoLGrid_clear_noinline (a);
	
	char line [65536];
	s32 x = 0;
	s32 y = 0;
	s32 run = 0;
	int finished = FALSE;
	while (!finished && fgets (line, sizeof (line), stdin))
	{
		if (line [0] == '#' || line [0] == 'x')
			continue;
		const char *c;
		for (c = line; *c && !finished; c++)
		{
			if (*c >= '0' && *c <= '9')
				run = 10 * run + (*c - '0');
			else if (*c == '$')
			{
				y += (run ? run : 1);
				x = 0;
				run = 0;
			}
			else if (*c == '!')
				finished = TRUE;
			else if (*c == 'b' || *c == '.')
			{
				x += (run ? run : 1);
				run = 0;
			}
			else if ((*c >= 'a' && *c <= 'z') || (*c >= 'A' && *c <= 'Z'))
			{
				s32 k;
				for (k = 0; k < (run ? run : 1); k++, x++)
					if (!GoLGrid_set_cell_on (a, x, y))
					{
						fprintf (stderr, "Pattern does not fit in the grid\n");
						return FALSE;
					}
				run = 0;
			}
		}
	}
	
	s64 gen;
	for (gen = 0; gen < gens; gen++)
	{
		GoLGrid_evolve (a, b);
		GoLGrid *t = a;
		a = b;
		b = t;
	}
	
	for (y = -(GRID_SIZE / 2); y < GRID_SIZE / 2; y++)
		for (x = -(GRID_SIZE / 2); x < GRID_SIZE / 2; x++)
			if (GoLGrid_get_cell (a, x, y))
				printf ("%d %d\n", x, y);
	
	Rect bb;
	GoLGrid_get_bounding_box (a, &bb);
	fprintf (stderr, "population %" PRIu64 ", bounding box %d %d %d %d\n", GoLGrid_get_population (a), bb.left_x, bb.top_y, bb.width, bb.height);
	return TRUE;
}

int main (int argc, char **argv)
{
	const char *rule_string = "B3/S23";
	const char *list_file_name = NULL;
	int rule_info = FALSE;
	int emit_c = FALSE;
	s64 evolve_gens = -1;
	int seeded = FALSE;
	s64 seed = 0;
	s64 thread_cnt = 1;

	s32 arg_ix;
	for (arg_ix = 1; arg_ix < argc; arg_ix++)
	{
		const char *a = argv [arg_ix];
		const char *v = (arg_ix + 1 < argc ? argv [arg_ix + 1] : NULL);
		s64 n = 0;
		int ok = TRUE;

		if (strcmp (a, "--evolve") == 0 && v)
		{
			if (!parse_s64_arg (v, 0, 1000000, &evolve_gens))
			{
				fprintf (stderr, "Bad value for --evolve: %s\n", v);
				return EXIT_FAILURE;
			}
			arg_ix++;
			continue;
		}
		if (strcmp (a, "--emit-c") == 0)
		{
			emit_c = TRUE;
			continue;
		}
		if (strcmp (a, "--rule-info") == 0)
		{
			rule_info = TRUE;
			continue;
		}
		if (strcmp (a, "-h") == 0 || strcmp (a, "--help") == 0)
		{
			print_usage ();
			return EXIT_SUCCESS;
		}
		if (a [0] != '-' || a [1] == '\0' || a [2] != '\0' || !v)
		{
			print_usage ();
			return EXIT_FAILURE;
		}

		switch (a [1])
		{
			case 'r': rule_string = v; break;
			case 'l': list_file_name = v; break;
			case 'n': ok = parse_s64_arg (v, 1, ((s64) 1) << 62, &used_soup_cnt); break;
			case 't': ok = parse_s64_arg (v, 1, ((s64) 1) << 62, &report_threshold); break;
			case 'k': ok = parse_s64_arg (v, 1, 1000000000, &max_key_patterns); break;
			case 'g': ok = parse_s64_arg (v, 0, 1000000000, &n); max_gens = n; break;
			case 'a': ok = parse_s64_arg (v, 1, MAX_STATE_A_POP_THRESHOLD, &n); state_a_pop_threshold = n; break;
			case 'b': ok = parse_s64_arg (v, 1, MAX_STATS_POP, &n); state_b_pop_threshold = n; break;
			case 'c': ok = parse_s64_arg (v, 1, MAX_STATS_POP, &n); state_c_unique_pop_threshold = n; break;
			case 's': ok = parse_s64_arg (v, 0, INT64_MAX, &seed); seeded = TRUE; break;
			case 'j': ok = parse_s64_arg (v, 1, 256, &thread_cnt); break;
			default: ok = FALSE; break;
		}

		if (!ok)
		{
			fprintf (stderr, "Bad value for %s: %s\n", a, v);
			return EXIT_FAILURE;
		}
		arg_ix++;
	}

	if (state_b_pop_threshold <= state_a_pop_threshold)
	{
		fprintf (stderr, "-b must be larger than -a\n");
		return EXIT_FAILURE;
	}

	#ifndef _OPENMP
		if (thread_cnt > 1)
		{
			fprintf (stderr, "Built without OpenMP; running with 1 thread\n");
			thread_cnt = 1;
		}
	#endif

	u8 table [512];
	if (!Rule_parse (rule_string, table))
		return EXIT_FAILURE;

	if (rule_info)
	{
		char map [90];
		Rule_map_string (table, map);
		Rule_init (&the_rule, table);
		printf ("%s\nBDD nodes: %d\nB0: %s\nIsotropic: %s\n", map, the_rule.node_cnt, table [0] ? "yes" : "no", Rule_table_is_isotropic (table) ? "yes" : "no");
		return EXIT_SUCCESS;
	}

	if (table [0])
	{
		fprintf (stderr, "B0 rules are not supported: the census grid has a fixed dead background\n");
		return EXIT_FAILURE;
	}


	if (!Rule_init (&the_rule, table))
		return EXIT_FAILURE;
	
	if (emit_c)
	{
		Rule_emit_c (&the_rule, stdout);
		return EXIT_SUCCESS;
	}
	
	#ifdef RULE_GENERATED
	{
		char map [90];
		Rule_map_string (table, map);
		if (strcmp (map, RULE_GENERATED_MAP) != 0)
		{
			fprintf (stderr, "This executable was built for the rule %s only; rebuild with make RULE=... or use the generic build\n", RULE_GENERATED_MAP);
			return EXIT_FAILURE;
		}
	}
	#endif


	if (evolve_gens >= 0)
		return evolve_test (evolve_gens) ? EXIT_SUCCESS : EXIT_FAILURE;
	
	if (!Rule_table_is_isotropic (table))
	{
		fprintf (stderr, "Non-isotropic rules are not supported: the census merges all 8 orientations of each pattern\n");
		return EXIT_FAILURE;
	}
	
	if (seeded)
	{
		u64 x = (u64) seed;
		u64 s0 = splitmix64 (&x);
		u64 s1 = splitmix64 (&x);
		random_u64_set_seed (s0, s1, FALSE);
	}

	if (!main_do ((s32) thread_cnt, list_file_name))
		return EXIT_FAILURE;

	return EXIT_SUCCESS;
}
