#define _GNU_SOURCE
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

char *USAGE = " \n\
Usage: %s OLD_DIR NEW_DIR DETAIL_OUT CLASS_COUNT_OUT CHANGED_READS_OUT \n\
                                                                       \n\
OLD_DIR should be the directory containing the data produced by the    \n\
old version of RUM, and NEW_DIR should contain the data produced by    \n\
the new version.                                                       \n\
                                                                       \n\
DETAIL_OUT will contain three columns: reads, old_alns, and            \n\
new_alns. Each row indicates the number of 'reads' that had 'old_alns' \n\
alignments in the old version and 'new_alns' alignments in the new     \n\
version.                                                               \n\
                                                                       \n\
CLASS_COUNT_OUT will contain three columns: reads, old_class,          \n\
new_class. A 'class' indicates the type of mapping: 'none', 'unique',  \n\
or 'non-unique'. Each row indicates that the number of reads that were \n\
in 'old_class' in the old version and 'new_class' in the new version.  \n\
                                                                       \n\
CHANGED_READS_OUT will contain a list of sequence ids that had a       \n\
different number of mappings on the old and new version.               \n\
";

// Largest read sequence number we expect to see. 
#define MAX_READS 100000000

// Largest number of alignments for a read that we want to capture.
#define MAX_ALNS 1000

// OLD_READS maps sequence number to the number of alignments in the
// "old" job, NEW_READS does the same for the "new" job.
int OLD_READS[MAX_READS];
int NEW_READS[MAX_READS];
int MATRIX[MAX_ALNS][MAX_ALNS];

const int DIR_NONE = 0;
const int DIR_FWD  = 1;
const int DIR_REV  = 2;

/**
 * Populate dest with a mapping from sequence number to the number of
 * alignments for the corresponding read, based on the contents of
 * filename, which should point to a RUM_Unique or RUM_NU file. max_id
 * will be adjusted to be the larger of the incoming max_id and the
 * largest sequence number read from the file.
 */
int load_reads_from_file(char *filename, int *dest, long *max_id) {
  printf("  Loading reads from %s\n", filename);

  FILE *f = fopen(filename, "r");
  if (!f) {
    perror(filename);
    return -1;
  }

  const int max_id_size = 20;

  char *line = NULL;
  size_t linecap = 0;
  long linelen = 0;

  int last_id  = 0;
  int last_dir = 0;

  char id_str[max_id_size];

  long count = 0;

  while ((linelen = getline(&line, &linecap, f)) > 0) {

    char *id_start_p = line + 4;
    int idlen = strcspn(id_start_p, "ab\t\n");
    
    int dir = DIR_NONE;
    
    switch (id_start_p[idlen]) {
      case 'a': dir = DIR_FWD; break;
      case 'b': dir = DIR_REV; break;
    }

    strncpy(id_str, id_start_p, idlen);
    id_str[idlen] = 0;
    long id = atol(id_str);
    
    if (id > MAX_READS) {
      fprintf(stderr, "ID too large: %ld\n", id);
    }
    if (id == last_id &&
        last_dir == DIR_FWD &&
        dir      == DIR_REV) {
    }
    else {
      dest[id]++;
    }

    last_id = id;
    last_dir = dir;

    if (id > *max_id) {
      *max_id = id;
    }
    
    if ((++count % 10000000) == 0) {
      printf("%10ld\n", count);
    }

  }

  fclose(f);
  printf("Max id is %ld\n", *max_id);
  return 0;
}

int load_reads_from_dir(char *dir, int *dest, long *max_id) {
  printf("Loading reads from %s\n", dir);

  char *suffixes[] = { "/RUM_Unique", 
                      "/RUM_NU" };

  int i;

  for (i = 0; i < 2; i++) {
    int len = strlen(dir) + strlen(suffixes[i]);
    char *filename = malloc(len);
    strcpy(filename, dir);
    strcat(filename, suffixes[i]);
    if (load_reads_from_file(filename, dest, max_id) < 0) {
      return -1;
    }
    free(filename);
  }
  
  return 0;
}

int main (int argc, char **argv) {

  if (argc < 6) {
    printf(USAGE, argv[0]);
    return 1;
  }

  char *old_dir          = argv[1];
  char *new_dir          = argv[2];
  char *out              = argv[3];
  char *changed_filename = argv[4];
  char *classes_filename = argv[5];

  long max_id = 0;

  printf("Initializing matrices\n");
  bzero(OLD_READS, sizeof(OLD_READS));
  bzero(NEW_READS, sizeof(NEW_READS));
  bzero(MATRIX, sizeof(MATRIX));

  if (load_reads_from_dir(old_dir, OLD_READS, &max_id) < 0) {
    printf("Exiting");
    return -1;
  }
  if (load_reads_from_dir(new_dir, NEW_READS, &max_id) < 0) {
    printf("Exiting");
    return -1;
  }

  FILE *changed_f = fopen(changed_filename, "w");
  if (!changed_f) { 
    perror(changed_filename); 
    return -1;
  };

  FILE *out_f = fopen(out, "w");  
  if (!out_f) { 
    perror(out); 
    return -1;
  };

  FILE *classes_f = fopen(classes_filename, "w");  
  if (!classes_f) { 
    perror(classes_filename); 
    return -1;
  };

  int aln_type_counts[3][3];
  bzero(aln_type_counts, sizeof(aln_type_counts));

  int i, j;
  fprintf(changed_f, "%s\t%s\t%s\n", "seqnum", "old_alns", "new_alns");

  for (i = 0; i < max_id; i++) {
    int old_reads = OLD_READS[i];
    int new_reads = NEW_READS[i];
    if (old_reads > MAX_ALNS ||
        new_reads > MAX_ALNS) {
      fprintf(stderr, 
              "Read %d has too many alignments (%d and %d), max is %d\n",
              i, old_reads, new_reads, MAX_ALNS);
    }
    else if (old_reads || new_reads) {
      MATRIX[old_reads][new_reads]++;
    }

    int old_aln_type 
      = old_reads == 0 ? 0 
      : old_reads == 1 ? 1 
      :                  2;
    int new_aln_type 
      = new_reads == 0 ? 0 
      : new_reads == 1 ? 1 
      :                  2;
    aln_type_counts[old_aln_type][new_aln_type]++;
    
    if (old_reads != new_reads) {
      fprintf(changed_f, "%d\t%d\t%d\n", i, old_reads, new_reads);
    }
  }

  // Print the class count summary file
  fprintf(classes_f, "%s\t%s\t%s\n", 
          "reads", "old_class", "new_class");
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      fprintf(classes_f, "%d\t%d\t%d\n",
              aln_type_counts[i][j], i, j);
    }
  }

  // Print the main counts output file
  fprintf(out_f, "%s\t%s\t%s\n", 
          "reads", "old_alns", "new_alns");
  for (i = 0; i < MAX_ALNS; i++) {
    for (j = 0; j < MAX_ALNS; j++) {
      int num_reads = MATRIX[i][j];
      if (num_reads) {
        fprintf(out_f, "%d\t%d\t%d\n", num_reads, i, j);
      }
    }
  }

  return 0;
}

