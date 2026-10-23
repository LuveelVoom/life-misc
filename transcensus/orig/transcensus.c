#include <stdlib.h>
#include <inttypes.h>
#include <memory.h>
#include <time.h>
#include <string.h>
#include <stdio.h>

#include "lib/lib.c"
#include "lib/rect.c"
#include "lib/celllist.c"
#include "lib/objects.c"
#include "lib/randomarray.c"
#include "lib/golgrid.c"
#include "lib/gridmisc.c"
#include "lib/hashtable.c"
#include "lib/golutils.c"
#include "lib/goltext.c"

#define GRID_SIZE 192
#define USED_SOUP_CNT ((s64) 10000000)
#define MAX_KEY_PATTERNS 2000000
#define MAX_MIN_PATTERNS 16
#define REPORT_THRESHOLD 200

#define STATE_A_POP_THRESHOLD 12
#define STATE_B_POP_THRESHOLD 20
#define STATE_C_UNIQUE_POP_THRESHOLD 4
#define MAX_STATS_POP 96
#define MAX_GENS 256

typedef struct
{
	s64 cnt;
	u64 hash;
	ObjCellList_Cell min_pattern_cell [STATE_A_POP_THRESHOLD];
	ObjCellList min_pattern;
} MinPattern;

typedef struct
{
	s64 cnt;
	s32 min_pattern_cnt;
	MinPattern min_pattern [MAX_MIN_PATTERNS];
} KeyPattern;

typedef enum {STATE_PRE_A, STATE_A, STATE_B, STATE_C} ObjectIdentificationState;

static const s64 min_value_multiplier [STATE_A_POP_THRESHOLD + 1] = {128, 128, 128, 128, 128, 96, 48, 24, 12, 7, 5, 4, 3};


static __not_inline GoLGrid *alloc_grid (void)
{
	return GoLUtils_alloc_std_grid (-(GRID_SIZE / 2), -(GRID_SIZE / 2), GRID_SIZE, GRID_SIZE);
}

static __not_inline int run_pattern (const GoLGrid *start_pattern, GoLGrid *key_pattern, GoLGrid *min_pattern)
{
	static int static_init = FALSE;
	static GoLGrid *p_m1;
	static GoLGrid *p_m0;
	
	if (!static_init)
	{
		static_init = TRUE;
		p_m1 = alloc_grid ();
		p_m0 = alloc_grid ();
	}
	
	GoLGrid_copy_noinline (start_pattern, p_m0);
	
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
		s32 pop = GoLGrid_get_population (p_m0);
		
		if (state == STATE_PRE_A)
			if (pop <= STATE_A_POP_THRESHOLD)
				state = STATE_A;
		
		if (state == STATE_A)
		{
			if (min_pop == -1 || pop < min_pop)
			{
				min_pop = pop;
				GoLGrid_copy (p_m0, min_pattern);
			}
			
			if (pop >= STATE_B_POP_THRESHOLD)
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
			if (pop <= STATE_A_POP_THRESHOLD)
				state = STATE_PRE_A;
			else
				if (pop < MAX_STATS_POP)
				{
					if (pop_stats [pop] == 0)
						unique_pops++;
					
					pop_stats [pop]++;
					
					if (unique_pops >= STATE_C_UNIQUE_POP_THRESHOLD)
					{
						GoLGrid_copy (p_m0, key_pattern);
						return TRUE;
					}
				}
		}
		
		if (gen >= MAX_GENS)
			return FALSE;
		
		GoLGrid *temp = p_m1;
		p_m1 = p_m0;
		p_m0 = temp;
		
		GoLGrid_evolve (p_m1, p_m0);
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

static __not_inline void store_min_pattern (const GoLGrid *min_pattern, u64 min_hash, KeyPattern *key_pattern)
{
	s32 min_ix;
	for (min_ix = 0; min_ix < key_pattern->min_pattern_cnt; min_ix++)
		if (key_pattern->min_pattern [min_ix].hash == min_hash)
		{
			key_pattern->min_pattern [min_ix].cnt++;
			return;
		}
	
	if (key_pattern->min_pattern_cnt >= MAX_MIN_PATTERNS)
		if (!remove_oldest_single (key_pattern))
			return;
	
	key_pattern->min_pattern [key_pattern->min_pattern_cnt].cnt = 1;
	key_pattern->min_pattern [key_pattern->min_pattern_cnt].hash = min_hash;
	
	ObjCellList_make_empty (&key_pattern->min_pattern [key_pattern->min_pattern_cnt].min_pattern, key_pattern->min_pattern [key_pattern->min_pattern_cnt].min_pattern_cell, STATE_A_POP_THRESHOLD);
	GoLGrid_to_obj_cell_list (min_pattern, &key_pattern->min_pattern [key_pattern->min_pattern_cnt].min_pattern);
	
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

static __not_inline int main_do (void)
{
	GoLGrid *soup = alloc_grid ();
	GoLGrid *key_pattern = alloc_grid ();
	GoLGrid *min_pattern = alloc_grid ();
	GoLGrid *canonical_key_pattern = alloc_grid ();
	GoLGrid *canonical_min_pattern = alloc_grid ();
	GoLGrid *canonical_t1 = alloc_grid ();
	GoLGrid *canonical_t2 = alloc_grid ();
	
	s32 key_pattern_cnt = 0;
	KeyPattern *key_pattern_store = malloc (MAX_KEY_PATTERNS * sizeof (KeyPattern));
	KeyPattern **sorted_key_pattern = malloc (MAX_KEY_PATTERNS * sizeof (KeyPattern *));
	
	RandomDataArray rda;
	RandomDataArray_create (&rda, GRID_SIZE * GRID_SIZE / 64);
	
	HashTable_u64 seen_key_patterns;
	HashTable_u64_create (&seen_key_patterns, 64, 0.75, 0.95);
	
	s64 settled_cnt = 0;
	s64 soup_ix = 0;
	while (TRUE)
	{
		if (settled_cnt >= USED_SOUP_CNT)
			break;
		
		if (soup_ix % 100000 == 0)
			fprintf (stderr, "Testing soup %" PRIi64 ", %" PRIi64 " used, %d unique key patterns\n", soup_ix, settled_cnt, key_pattern_cnt);
		
		soup_ix++;
		
		GoLGrid_clear_noinline (soup);
		GoLGrid_or_8_by_8_block_noinline (soup, 0, 0, random_u64 ());
		
		if (run_pattern (soup, key_pattern, min_pattern))
		{
			settled_cnt++;
			
			u64 key_hash = GoLUtils_make_canonical (key_pattern, canonical_key_pattern, &rda, canonical_t1, canonical_t2);
			u64 min_hash = GoLUtils_make_canonical (min_pattern, canonical_min_pattern, &rda, canonical_t1, canonical_t2);
			
			u64 key_pattern_ix;
			if (HashTable_u64_get_data (&seen_key_patterns, key_hash, &key_pattern_ix))
				key_pattern_store [key_pattern_ix].cnt++;
			else
			{
				if (key_pattern_cnt >= MAX_KEY_PATTERNS)
				{
					fprintf (stderr, "Key pattern store overflow\n");
					return FALSE;
				}
				
				key_pattern_ix = key_pattern_cnt;
				key_pattern_cnt++;
				
				HashTable_u64_store (&seen_key_patterns, key_hash, key_pattern_ix, TRUE, NULL);
				
				key_pattern_store [key_pattern_ix].cnt = 1;
				key_pattern_store [key_pattern_ix].min_pattern_cnt = 0;
			}
			
			store_min_pattern (canonical_min_pattern, min_hash, &key_pattern_store [key_pattern_ix]);
		}
	}
	
	create_sorted_key_patterns (key_pattern_cnt, key_pattern_store, sorted_key_pattern);
	
	GoLGrid *out_grid_p = GoLUtils_alloc_std_grid (0, 0, 2048, 16384);
	GoLGrid *out_grid_t = GoLUtils_alloc_std_grid (0, 0, 2048, 16384);
	
	s32 x_pos = 64;
	s32 y_pos = 64;
	
	s32 store_ix;
	s32 min_ix;
	for (store_ix = 0; store_ix < key_pattern_cnt; store_ix++)
		if (sorted_key_pattern [store_ix]->cnt >= REPORT_THRESHOLD)
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
			
			char buf [16];
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
	
	fprintf (stderr, "Used store positions = %d\n", key_pattern_cnt);
	fprintf (stderr, "Settled: %" PRIi64 " of %" PRIi64 "\n", settled_cnt, soup_ix);
	
	return TRUE;
}

int main (void)
{
	if (!verify_cpu_type ())
		return EXIT_FAILURE;
	
	if (!main_do ())
		return EXIT_FAILURE;
	
	return EXIT_SUCCESS;
}
