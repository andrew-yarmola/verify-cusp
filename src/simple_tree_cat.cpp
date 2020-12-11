#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <string.h>
#include <stdlib.h>
#include <dirent.h>
#include <vector>
#include <string>

bool g_recursive = false;
bool g_verbose = false;
bool g_start_is_root = false;
char* g_tree_location;

FILE* open_box(char* boxcode)
{
  FILE* fp;
  char fileboxcode[10000];
  char file_name[10000];
  strcpy(fileboxcode, boxcode);
  // Open the root file if empty
  if (strcmp(fileboxcode, "") == 0)
    strcpy(fileboxcode, "root");
  sprintf(file_name, "%s/%s.out", g_tree_location, fileboxcode);
  struct stat sb;
  if (0 == stat(file_name, &sb)) {
    if (g_verbose) fprintf(stderr, "opening %s\n", file_name);
    fp = fopen(file_name, "r");
    return fp;
  }
  // Look for a gzipped file
  sprintf(file_name, "%s/%s.out.tar.gz", g_tree_location, fileboxcode);
  if (0 == stat(file_name, &sb)) {
    char command_buf[10000];
    if (g_verbose) fprintf(stderr, "opening %s\n", file_name);
    sprintf(command_buf, "tar -xOzf %s", file_name);
    fp = popen(command_buf, "r");
    return fp;
  }
  return 0;
}

bool process_tree(FILE* fp, bool print_tree, char* boxcode)
{
  bool success = true;
  int boxdepth = strlen(boxcode);
  char buf[10000];
  int depth = 0;
  while (fgets(buf, sizeof(buf), fp)) {
    // Open HOLE file if exists. If no printing is set, we don't need to traverse HOLEs
    if (buf[0] == 'H') {
      if (g_recursive && depth > 0) {
        FILE* fph = open_box(boxcode);
        if (fph) {
          success = process_tree(fph, print_tree, boxcode);
          fclose(fph);
          if (!success) {
            fprintf(stderr, "Fatal! Hole file %s incomplete\n", boxcode);
            // The tree is incomplete, so we rename the boxfile to mark as incomplete
            break;
          }
        }
      }
    } else if (print_tree) {
      fprintf(stdout, "%s", buf); 
    }
    if (buf[0] == 'X') {
      boxcode[boxdepth + depth] = '0'; // Descend via left branch
      ++depth;
      boxcode[boxdepth + depth] = '\0';
    } else {
      // Go up as many nodes as necessary
      for (; depth > 0 && boxcode[boxdepth + depth-1] == '1'; --depth) {}
      if (depth > 0) {
        boxcode[boxdepth + depth-1] = '1'; // Jump from left to right node
        boxcode[boxdepth + depth] = '\0'; // Truncate to keep box current
      } else {
        boxcode[boxdepth] = '\0'; // Truncate to keep box current
        return true;
      }
    }
  }
  return false; 
}

int main(int argc, char** argv)
{
  bool print_tree = true;

  if (argc > 1 && strcmp(argv[1], "-s") == 0) {
    print_tree = false;
    ++argv;
    --argc;
  }

  if (argc > 1 && strcmp(argv[1], "-v") == 0) {
    g_verbose = true;
    ++argv;
    --argc;
  }

  if (argc > 1 && strcmp(argv[1], "-r") == 0) {
    g_recursive = true;
    ++argv;
    --argc;
  }

  if (argc != 3) {
    fprintf(stderr, "Usage: simple_tree_cat [-s] [-v] [-r] tree_location boxcode\n");
    fprintf(stderr, "s : silent, don't print trees or holes\n");
    fprintf(stderr, "v : verbose\n");
    fprintf(stderr, "r : recur over all subtree files to prince full subtree of boxcode\n");
    exit(1);
  }

  // The fullboxcode parameter can specify the file_name and sequetial boxcode
  // A boxcode is just a sequence of zeros and ones giving a posiiton in a binary tree depth-first traversal
  // The tree file will also be in pre-order depth-first
  char fullboxcode[10000];
  char fileboxcode[10000];
  strncpy(fullboxcode, argv[2], 10000);
  strncpy(fileboxcode, argv[2], 10000);
  g_tree_location = argv[1];

  int file_box_len = strlen(fullboxcode);

  if ((file_box_len == 0) || (strncmp(fullboxcode, "root", file_box_len) == 0)) {
    g_start_is_root = true;
  }

  // See if a file with the tree for a prefix of the box exists
  FILE* fp = 0;
  while (file_box_len >= 0) {
    fileboxcode[file_box_len] = '\0';
    fp = open_box(fileboxcode);
    if (fp) { break; }
    --file_box_len;
  }

  if (!fp) { 
    fprintf(stderr, "Fatal! Failed to open boxfile\n");
    exit(2);
  }

  char * boxcode_const = (char *)calloc(10000, sizeof(char));
  strncpy(boxcode_const, fullboxcode+file_box_len, 10000);
  char * boxcode = boxcode_const;

  char buf[10000];
  // If the boxcode is still not empty, we traverse down the tree and print only
  // once we get to the proper node. We terminate early if the node does not exist
  while (*boxcode && fgets(buf, sizeof(buf), fp)) {
    if (buf[0] != 'X') { // If not a splitting, print the test failed by the truncated box
      *boxcode = '\0';
      //fprintf(stderr, "terminal box = %s%s\n", fileboxcode, boxcode_const);
      if (print_tree) fputs(buf, stdout);
      free(boxcode_const);
      fclose(fp);
      exit(0);
    }
    if (*boxcode == '1') { // Actually have to process the tree if we go right at any point in the boxcode
      bool success = process_tree(fp, false, boxcode);
      if (!success) {
        fprintf(stderr, "Fatal! Incomplete tree\n");
        exit(2); // Incomplete tree or boxcode not found
      }
    }
    ++boxcode; // Keeps going left in the tree as *boxcode == 0
  }

  free(boxcode_const);

  bool success = process_tree(fp, print_tree, fullboxcode);
  fclose(fp);

  if (!success) {
    // The tree is incomplete, so we rename the boxfile to mark as incomplete
    // TODO: Not sure if treecat should have the power to rename files
    fprintf(stderr, "Fatal! Incomplete tree\n");
    exit(2);
  }
  exit(0); 
}
