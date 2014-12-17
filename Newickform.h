#ifndef __NEWICKFORM_H__
#define __NEWICKFORM_H__

// NOTE: this code was downloaded from http://yuweibioinfo.blogspot.com/2008/10/newick-tree-parser-in-c-make-use-of.html

#include <list>
#include <set>
#include <string>
#include <vector>
#include <map>


typedef struct newick_child
{
	struct newick_node *node;
	struct newick_child *next;
} newick_child;

typedef struct newick_node
{
	char *taxon;
	char *seq;
	char* uniq;
	float dist;
	float height;
	int childNum;
	struct newick_child *child;
	struct newick_node *parent;
	std::vector<int> descendants;
} newick_node;

#ifdef __NEWICKFORM_C__
newick_node* parseTree(char *str);
void printTree(newick_node *root);
void getLeaves(newick_node *root, std::list<newick_node*>& L);
void removeDuplicates(std::vector<int>& vec);
void getCommon(char * pcTreeStr, std::string **);
int NUM;

#else
extern newick_node* parseTree(char *str);
extern void printTree(newick_node *root);
extern void getLeaves(newick_node *root);
extern void getLeaves(newick_node *root, std::list<newick_node*>& L);
extern void removeDuplicates(std::vector<int>& vec);
extern void getCommon(char * pcTreeStr, std::string **);
extern int NUM;

#endif

#endif

