// NOTE: this code was downloaded from http://yuweibioinfo.blogspot.com/2008/10/newick-tree-parser-in-c-make-use-of.html

#define __NEWICKFORM_C__

#include "seqUtil.h"
#include "Newickform.h"
#include <vector>
#include <list>
#include <iostream>
#include <algorithm>
#include <sstream>

using namespace std;

void getCommon(char * pcTreeStr, std::string ** common) {

//	cout << "received " << pcTreeStr << endl;
	
        newick_node *root;

	// Parse tree string
	root = parseTree(pcTreeStr);
//	printTree(root);
//	printf("\n");
	std::list<newick_node*> L;
	getLeaves(root, L);
	int N = L.size();
//	cout << "ELEMENTS: " << N << endl;

        for (int i = 0; i<N; i++) {
		for (int j = 0; j<N; j++) {
                	common[i][j] = "";
		}
        }


//        for(list<newick_node*>::iterator i=L.begin(); i != L.end(); ++i) {
//                if ((*i)->parent != NULL) cout << "PARENT: " << string((*i)->parent->taxon) << endl;
//        }

//        for(list<newick_node*>::iterator i=L.begin(); i != L.end(); ++i) {
//                cout << "SES " << string((*i)->taxon) << endl;
//        }
	
	int counter = 0;
        while (L.size() > 0) { //for each node in the todo list, climb up the tree and build the descendant sets.
            counter++;
            newick_node* node = L.front();
	    L.pop_front();

	    // create unique ID
            stringstream sst;
            sst << node->taxon << "_" << 0+node->height;
//	    node->uniq = sst.str().c_str();
	    node->uniq = new char[sst.str().size() + 1]; 
	    strcpy(node->uniq, sst.str().c_str());

	    string name = string(node->taxon);
	    string parent;
 	    if (node->parent != NULL) {
		parent = string(node->parent->taxon);
         	node->parent->height = node->dist+node->height;
                L.push_back(node->parent);
		std::vector<int> v(node->parent->descendants.size() + node->descendants.size());
		merge(node->parent->descendants.begin(), node->parent->descendants.end(), node->descendants.begin(), node->descendants.end(), v.begin());
		removeDuplicates(v);
		node->parent->descendants = v;
 	    }
	    else {
		parent = "ROOT";
	    }
//            cout << "node: " << name << " parent: " << parent << " distance: " << node->dist << " height " << node->height << " uniq " << node->uniq << " descendants: ";
//	    for (int i=0; i<node->descendants.size(); i++) cout << node->descendants[i] << " ";
//		cout << endl;
	}

        L.clear();
        counter = 0;
        L.push_back(root);
        while (L.size() > 0) {
            counter++;
            newick_node* node = L.front();
            L.pop_front();

	    for (int c1 = 1; c1 <= node->childNum; c1++) {
		newick_child* child1 = node->child;
		for (int count = 2; count <= c1; count++) { child1 = child1->next; }
//		cout << child1->node->uniq << endl;
		newick_node *n1 = child1->node;
                L.push_back(n1);
//		cout << "N1: " << n1->uniq << endl;
		for (int c2 = c1+1; c2 <= node->childNum; c2++) {
	                newick_child* child2 = node->child;
	                for (int count = 2; count <= c2; count++) { child2 = child2->next; }
	                newick_node *n2 = child2->node;
//	                cout <<	"N1: " << n1->uniq << " N2: " << n2->uniq << endl;
			for (int i = 0; i < n1->descendants.size(); i++) {
				for (int j = 0; j < n2->descendants.size(); j++) {
					int d1 = n1->descendants[i];
                                        int d2 = n2->descendants[j];
//					cout << d1 << " vs " << d2 << " common: " << node->uniq << endl;
					if (d2 > d1) common[d1-1][d2-1] = node->uniq;
                                        else common[d2-1][d1-1]	= node->uniq;
				}
			}
		}
	    }
	}

	// Free occupied memory
//	seqFree(pcOutputFile);
	seqFree(pcTreeStr);

	// End memory management procedure and free all allocated space
	seqFreeAll();

	NUM = N;

//	return common;
}



newick_node* parseTree(char *str)
{
	newick_node *node;
	newick_child *child;
	char *pcCurrent;
	char *pcStart;
	char *pcColon = NULL;
	char cTemp;
	int iCount;

//	if (sDEBUG == 1)
//	{
//		printf("%s\n", str);
//	}
	pcStart = str;

	if (*pcStart != '(')
	{
		// Leaf node. Separate taxon name from distance. If distance not exist then take care of taxon name only
		pcCurrent = str;
		while (*pcCurrent != '\0')
		{
			if (*pcCurrent == ':')
			{
				pcColon = pcCurrent;
			}
			pcCurrent++;
		}
		node = (newick_node*)seqMalloc(sizeof(newick_node));
		if (pcColon == NULL)
		{
			// Taxon only
			node->taxon = (char*)seqMalloc(strlen(pcStart) + 1);
			memcpy(node->taxon, pcStart, strlen(pcStart));
		}
		else
		{
			// Taxon
			*pcColon = '\0';
			node->taxon = (char*)seqMalloc(strlen(pcStart) + 1);
			memcpy(node->taxon, pcStart, strlen(pcStart));
			*pcColon = ':';
			// Distance
			pcColon++;
			node->dist = (float)atof(pcColon);
		}
		node->childNum = 0;
	}
	else
	{
		// Create node
		node = (newick_node*)seqMalloc(sizeof(newick_node));
		child = NULL;
		// Search for all child nodes
		// Find all ',' until corresponding ')' is encountered
		iCount = 0;
		pcStart++;
		pcCurrent = pcStart;
		while (iCount >= 0)
		{
			switch (*pcCurrent)
			{
				case '(':
					// Find corresponding ')' by counting
					pcStart = pcCurrent;
					pcCurrent++;
					iCount++;
					while (iCount > 0)
					{
						if (*pcCurrent == '(')
						{
							iCount++;
						}
						else if (*pcCurrent == ')')
						{
							iCount--;
						}
						pcCurrent++;
					}
					while (*pcCurrent != ',' && *pcCurrent != ')')
					{
						pcCurrent++;
					}
 					cTemp = *pcCurrent;
					*pcCurrent = '\0';
					// Create a child node
					if (child == NULL)
					{
						node->child = (newick_child*)seqMalloc(sizeof(newick_child));
						node->childNum = 1;
						child = node->child;
					}
					else
					{
						child->next = (newick_child*)seqMalloc(sizeof(newick_child));
						node->childNum++;
						child = child->next;
					}
					child->node = parseTree(pcStart);
					*pcCurrent = cTemp;
					if (*pcCurrent != ')')
					{
						pcCurrent++;
					}
				break;

				case ')':
					// End of tihs tree. Go to next part to retrieve distance
					iCount--;
				break;

				case ',':
					// Impossible separation since according to the algorithm, this symbol will never encountered.
					// Currently don't handle this and don't create any node
				break;

				default:
					// leaf node encountered
					pcStart = pcCurrent;
					while (*pcCurrent != ',' && *pcCurrent != ')')
					{
						pcCurrent++;
					}
					cTemp = *pcCurrent;
					*pcCurrent = '\0';
					// Create a child node
					if (child == NULL)
					{
						node->child = (newick_child*)seqMalloc(sizeof(newick_child));
						node->childNum = 1;
						child = node->child;
					}
					else
					{
						child->next = (newick_child*)seqMalloc(sizeof(newick_child));
						node->childNum++;
						child = child->next;
					}
					child->node = parseTree(pcStart);
					*pcCurrent = cTemp;
					if (*pcCurrent != ')')
					{
						pcCurrent++;
					}
				break;
			}
		}

		// If start at ':', then the internal node has no name.
		pcCurrent++;
		if (*pcCurrent == ':')
		{
			pcStart = pcCurrent + 1;
			while (*pcCurrent != '\0' && *pcCurrent != ';')
			{
				pcCurrent++;
			}
			cTemp = *pcCurrent;
			*pcCurrent = '\0';
			node->dist = (float)atof(pcStart);
			*pcCurrent = cTemp;
		}
		else if (*pcCurrent != ';' && *pcCurrent != '\0')
		{
			// Find ':' to retrieve distance, if any.
			// At this time *pcCurrent should equal to ')'
			pcStart = pcCurrent;
			while (*pcCurrent != ':' && *pcCurrent != ';')
			{
				pcCurrent++;
			}
			cTemp = *pcCurrent;
			*pcCurrent = '\0';
			node->taxon = (char*)seqMalloc(strlen(pcStart) + 1);
			memcpy(node->taxon, pcStart, strlen(pcStart));
			*pcCurrent = cTemp;
			pcCurrent++;
			pcStart = pcCurrent;
			while (*pcCurrent != '\0' && *pcCurrent != ';')
			{
				pcCurrent++;
			}
			cTemp = *pcCurrent;
			*pcCurrent = '\0';
			node->dist = (float)atof(pcStart);
			*pcCurrent = cTemp;
		}
	}

	return node;
}

void getLeaves(newick_node *root, std::list<newick_node*>& L) {
        newick_child *child;
        if (root->childNum == 0)
        {
//                printf("%s\n", root->taxon, root->dist);
		L.push_back(root);
		int d = atoi(root->taxon);
		root->descendants.push_back(d);
		removeDuplicates(root->descendants);
        }
        else
        {
                child = root->child;
                while (child != NULL)
                {
			child->node->parent = root;
//			cout << child->node->parent << endl;
                        getLeaves(child->node, L);
                        child = child->next;
                }
        }
}

    void removeDuplicates(std::vector<int>& vec)
    {
        std::sort(vec.begin(), vec.end());
        vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
    }

void printTree(newick_node *root)
{
	newick_child *child;
	if (root->childNum == 0)
	{
		printf("%s:%0.6f", root->taxon, root->dist);
	}
	else
	{
		child = root->child;
		printf("(");
		while (child != NULL)
		{
			printTree(child->node);
			if (child->next != NULL)
			{
				printf(",");
			}
			child = child->next;
		}
		if (root->taxon != NULL)
		{
			printf(")%s:%0.6f", root->taxon, root->dist);
		}
		else
		{
			printf("):%0.6f", root->dist);
		}
	}
}
