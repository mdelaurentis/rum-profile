#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define MAX_READS 100000000
#define MAX_ALNS 1000

int OLD_READS[MAX_READS];
int NEW_READS[MAX_READS];
int MATRIX[MAX_ALNS][MAX_ALNS];

const int DIR_NONE = 0;
const int DIR_FWD = 1;
const int DIR_REV = 2;

int load_reads_from_file(char *filename, int *dest, long *max_id) {
  printf("  Loading reads from %s\n", filename);

  FILE *f = fopen(filename, "r");

  const int max_id_size = 20;

  char *line = NULL;
  long linecap = 0;
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
      printf("ID too large: %ld\n", id);
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
      printf("%10ld: %ld, %ld, %s\n", count, id, linelen, line);
    }

  }

  fclose(f);
  printf("Max id is %ld\n", *max_id);
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
    load_reads_from_file(filename, dest, max_id);
    free(filename);
  }
  
  return 0;
}

int main (int argc, char **argv) {
  if (argc < 4) {
    printf("Usage: %s OLD_DIR NEW_DIR OUT_FILE\n", argv[0]);
    return 1;
  }

  char *old_dir = argv[1];
  char *new_dir = argv[2];
  char *out     = argv[3];

  long max_id = 0;

  bzero(OLD_READS, sizeof(OLD_READS));
  bzero(NEW_READS, sizeof(NEW_READS));
  bzero(MATRIX, sizeof(MATRIX));
  load_reads_from_dir(old_dir, OLD_READS, &max_id);
  load_reads_from_dir(new_dir, NEW_READS, &max_id);

  printf("Max max is %ld\n", max_id);

  int max_alns = 0;

  printf("Max alns is %d\n", max_alns);

  int i, j;
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
  }

  FILE *out_f = fopen(out, "w");

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
