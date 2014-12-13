#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "kdtree.h"

#if defined(WIN32) || defined(__WIN32__)
#include <malloc.h>
#endif

#ifdef USE_LIST_NODE_ALLOCATOR

#ifndef NO_PTHREADS
#include <pthread.h>
#else

#ifndef I_WANT_THREAD_BUGS
#error "You are compiling with the fast list node allocator, with pthreads disabled! This WILL break if used from multiple threads."
#endif	/* I want thread bugs */

#endif	/* pthread support */
#endif	/* use list node allocator */

//using namespace std;

struct kdhyperrect {
	int dim;
	double *min, *max;              /* minimum/maximum coords */
};

struct kdnode {
	double *pos;
	int dir;
	int *data;
	struct kdnode *left, *right;	/* negative/positive side */
};

struct res_node {
	struct kdnode *item;
	double dist_sq;
	struct res_node *next;
};

struct kdtree {
	int dim;
	struct kdnode *root;
	struct kdhyperrect *rect;
	void (*destr)(void*);
};

struct kdres {
	struct kdtree *tree;
	struct res_node *rlist, *riter;
	int size;
};


#define SIZE 100

struct rheap {
    struct res_node* node;
    //double* heaparray;
    int capacity;
    int size;
};

#define SQ(x)			((x) * (x))


static void clear_rec(struct kdnode *node, void (*destr)(void*));
static int insert_rec(struct kdnode **node, const double *pos, int data, int dir, int dim);
static int rlist_insert(struct res_node *list, struct kdnode *item, double dist_sq);
static void clear_results(struct kdres *set);

static struct kdhyperrect* hyperrect_create(int dim, const double *min, const double *max);
static void hyperrect_free(struct kdhyperrect *rect);
static struct kdhyperrect* hyperrect_duplicate(const struct kdhyperrect *rect);
static void hyperrect_extend(struct kdhyperrect *rect, const double *pos);
static double hyperrect_dist_sq(struct kdhyperrect *rect, const double *pos);

#ifdef USE_LIST_NODE_ALLOCATOR
static struct res_node *alloc_resnode(void);
static void free_resnode(struct res_node*);
#else
#define alloc_resnode()		malloc(sizeof(struct res_node))
#define free_resnode(n)		free(n)
#endif



struct kdtree *kd_create(int k)
{
	struct kdtree *tree;

	if(!(tree = malloc(sizeof *tree))) {
		return 0;
	}

	tree->dim = k;
	tree->root = 0;
	tree->destr = 0;
	tree->rect = 0;

	return tree;
}

void kd_free(struct kdtree *tree)
{
	if(tree) {
		kd_clear(tree);
		free(tree);
	}
}

static void clear_rec(struct kdnode *node, void (*destr)(void*))
{
	if(!node) return;

	clear_rec(node->left, destr);
	clear_rec(node->right, destr);
	
	if(destr) {
		destr(node->data);
	}
	free(node->pos);
	free(node);
}

void kd_clear(struct kdtree *tree)
{
	clear_rec(tree->root, tree->destr);
	tree->root = 0;

	if (tree->rect) {
		hyperrect_free(tree->rect);
		tree->rect = 0;
	}
}

void kd_data_destructor(struct kdtree *tree, void (*destr)(void*))
{
	tree->destr = destr;
}


static int insert_rec(struct kdnode **nptr, const double *pos, int data, int dir, int dim)
{
	int new_dir;
	struct kdnode *node;

	if(!*nptr) {
		if(!(node = malloc(sizeof *node))) {
			return -1;
		}
		if(!(node->pos = malloc(dim * sizeof *node->pos))) {
			free(node);
			return -1;
		}
		memcpy(node->pos, pos, dim * sizeof *node->pos);
		node->data = data;
		node->dir = dir;
		node->left = node->right = 0;
		*nptr = node;
		return 0;
	}

	node = *nptr;
	new_dir = (node->dir + 1) % dim;
	if(pos[node->dir] < node->pos[node->dir]) {
		return insert_rec(&(*nptr)->left, pos, data, new_dir, dim);
	}
	return insert_rec(&(*nptr)->right, pos, data, new_dir, dim);
}

int kd_insert(struct kdtree *tree, const double *pos, int data)
{
	if (insert_rec(&tree->root, pos, data, 0, tree->dim)) {
		return -1;
	}

	if (tree->rect == 0) {
		tree->rect = hyperrect_create(tree->dim, pos, pos);
	} else {
		hyperrect_extend(tree->rect, pos);
	}

	return 0;
}
bool sortNode(struct kdnode* i,struct kdnode* j, int dim){
    printf("sorting\n");
    return i->pos[dim] < j->pos[dim];
}
static int insert_rec_all(struct kdnode *nptr, char* side, struct kdnode* nodes , int start, int end, int dim)
{
    if(!start == end){
        if ((end - start) == 1) {
            if (strcmp(side, "left")==0) {
                nptr->left=&nodes[start];
            }
            else if (strcmp(side, "right")==0) {
                nptr->right=&nodes[start];
            }
            return 0;
        }
        struct kdnode* new_node = &nodes[start];
        qsort_r(new_node, end-start,sizeof(nodes[0]),sortNode,dim);
        int median = start + (end-start)/2;
        if (strcmp(side, "left")==0) {
            nptr->left=&nodes[median];
            insert_rec_all(nptr->left, "left", nodes, start, median, dim+1);
            insert_rec_all(nptr->right, "right", nodes, median+1, end, dim+1);
        }
        else if (strcmp(side, "right")==0) {
           nptr->right=&nodes[median];
            insert_rec_all(nptr->left, "left", nodes, start, median, dim+1);
            insert_rec_all(nptr->right, "right", nodes, median+1, end, dim+1);
        }
    }
    return 0;
}



int kd_insert_all(struct kdtree *tree, const double **pos_all, int *data, int n, int D)
{   struct kdnode *nodes[n];
    int j;
    //printf("for\n");
    for (j=0; j<n; j++) {
        struct kdnode *node;
        //printf("for\n");
        if(!(node = malloc(sizeof *node))) {
            return -1;
        }
        if(!(node->pos = malloc(D * sizeof *node->pos))) {
            free(node);
            return -1;
        }
        //printf("for done\n");
        memcpy(node->pos, &pos_all[j], D * sizeof *node->pos);
        node->data = data[j];
        //node->dir = dir;
        node->left = node->right = 0;
        nodes[j] = node;
        free(node);
    }
    printf("for done\n");
    qsort_r(nodes, n, sizeof(nodes[0]),sortNode,0);
    printf("done sorting\n");
    tree->root = nodes[n/2];
    if (insert_rec_all(&tree->root,"left",&nodes, 0, n/2,1) && insert_rec_all(&tree->root, "right",&nodes,n/2+1 ,n, 1)) {
        return -1;
    }
    
    return 0;
}


static int find_nearest(struct kdnode *node, const double *pos, double range, struct res_node *list, int ordered, int dim)
{
	double dist_sq, dx;
	int i, ret, added_res = 0;

	if(!node) return 0;

	dist_sq = 0;
	for(i=0; i<dim; i++) {
		dist_sq += SQ(node->pos[i] - pos[i]);
	}
	if(dist_sq <= SQ(range)) {
		if(rlist_insert(list, node, ordered ? dist_sq : -1.0) == -1) {
			return -1;
		}
		added_res = 1;
	}

	dx = pos[node->dir] - node->pos[node->dir];

	ret = find_nearest(dx <= 0.0 ? node->left : node->right, pos, range, list, ordered, dim);
	if(ret >= 0 && fabs(dx) < range) {
		added_res += ret;
		ret = find_nearest(dx <= 0.0 ? node->right : node->left, pos, range, list, ordered, dim);
	}
	if(ret == -1) {
		return -1;
	}
	added_res += ret;

	return added_res;
}
struct rheap* initHeap() {
    struct rheap* h;
    h = (struct rheap*)(malloc(sizeof(struct rheap)));
    h->capacity = SIZE;
    //h->heaparray = (double*)malloc(sizeof(double)*(SIZE+1));
    h->size = 0;
    h->node=(struct res_node*)malloc(sizeof(struct res_node)*(SIZE+1));
    return h;
}
void swap(struct rheap *h, int index1, int index2) {
    struct res_node temp = h->node[index1];
    h->node[index1] = h->node[index2];
    h->node[index2] = temp;
}
void percolateUp(struct rheap *h, int index) {
    if (index > 1) {
        if (h->node[index/2].dist_sq < h->node[index].dist_sq) {
            swap(h, index, index/2);
            percolateUp(h, index/2);
        }
    }
}
int maximum(double a, int indexa, double b, int indexb) {
    if (a > b)
    return indexa;
    else
    return indexb;
}
void percolateDown(struct rheap *h, int index) {
    int max;
    if ((2*index+1) <= h->size) {
        max = maximum(h->node[2*index].dist_sq, 2*index, h->node[2*index+1].dist_sq, 2*index+1);
        if (h->node[index].dist_sq > h->node[max].dist_sq) {
            swap(h, index, max);
            percolateDown(h, max);
        }
    }
    else if (h->size == 2*index) {
        if (h->node[index].dist_sq > h->node[2*index].dist_sq)
        swap(h, index, 2*index);
    }
}

int rheap_insert(struct rheap *h, struct knode *node,double value) {
    //printf("inserting ");
    struct res_node* tempnode;
    struct res_node* thrownode;
    int i;
    if (h->size == h->capacity) {
        h->capacity *= 2;
        tempnode=(struct res_node*)malloc(sizeof(struct res_node)*h->capacity+1);
        for (i=0; i<h->capacity; i++){
            tempnode[i] = h->node[i];
        }
        thrownode= h->node;
        h->node = thrownode;
        free(thrownode);
    }
    h->size++;
    h->node[h->size].item=node;
    h->node[h->size].dist_sq=value;
    percolateUp(h, h->size);
    return 1;
}

struct res_node* rheap_remove_max(struct rheap *h) {
    struct res_node* retval;
    if (h->size > 0) {
        retval = &h->node[1];
        h->node[1] = h->node[h->size];
        h->size--;
        percolateDown(h, 1);
        return &retval;
    }
}
void heapify(struct rheap *h) {
    int i;
    for (i=h->size/2; i>0; i--)
    percolateDown(h, i);
    
}

struct rheap * initHeapfromArray(struct res_node* values, int length) {
    int i;
    struct rheap* h;
    h = (struct rheap*)(malloc(sizeof(struct rheap)));
    h->node = (struct res_node*)malloc(sizeof(struct res_node)*(length+1));
    for (i=0; i<length; i++)
    h->node[i] = values[i];
    h->size = length;
    heapify(h);
    return h;
}

void rheap_heapsort(struct res_node* values[], int length) {
    struct rheap *h;
    int i;
    h =  initHeapfromArray(values, length);
    length = h->size;
    for (i=0; i<length; i++) {
        values[i] = rheap_remove_max(h);
    }
}

struct res_node* rheap_get_max(struct rheap *h){
    return &h->node[1];
}
//#if 0
static int find_nearest_n(struct kdnode *node, const double *pos, double range, int num, struct rheap *heap, int dim)
{
	double dist_sq, dx;
	int i, ret, added_res = 0;
    double range_sq = SQ(range);
	if(!node) return 0;
    
	dist_sq = 0;
	for(i=0; i<dim; i++) {
		dist_sq += SQ(node->pos[i] - pos[i]);
	}
	if(dist_sq <= range_sq) {
		if(heap->size >= num) {
			/* get furthest element */
			struct res_node *maxelem = rheap_get_max(heap);

			/* and check if the new one is closer than that */
			if(maxelem->dist_sq > dist_sq) {
				rheap_remove_max(heap);

				if(rheap_insert(heap, node, dist_sq) == -1) {
					return -1;
				}
				added_res = 1;

				range_sq = dist_sq;
			}
		} else {
			if(rheap_insert(heap, node, dist_sq) == -1) {
				return -1;
			}
			added_res = 1;
		}
	}


	/* find signed distance from the splitting plane */
	dx = pos[node->dir] - node->pos[node->dir];

	ret = find_nearest_n(dx <= 0.0 ? node->left : node->right, pos, range, num, heap, dim);
	if(ret >= 0 && fabs(dx) < range) {
		added_res += ret;
		ret = find_nearest_n(dx <= 0.0 ? node->right : node->left, pos, range, num, heap, dim);
	}

}
//#endif

static void kd_nearest_i(struct kdnode *node, const double *pos, struct kdnode **result, double *result_dist_sq, struct kdhyperrect* rect)
{
	int dir = node->dir;
	int i;
	double dummy, dist_sq;
	struct kdnode *nearer_subtree, *farther_subtree;
	double *nearer_hyperrect_coord, *farther_hyperrect_coord;

	/* Decide whether to go left or right in the tree */
	dummy = pos[dir] - node->pos[dir];
	if (dummy <= 0) {
		nearer_subtree = node->left;
		farther_subtree = node->right;
		nearer_hyperrect_coord = rect->max + dir;
		farther_hyperrect_coord = rect->min + dir;
	} else {
		nearer_subtree = node->right;
		farther_subtree = node->left;
		nearer_hyperrect_coord = rect->min + dir;
		farther_hyperrect_coord = rect->max + dir;
	}

	if (nearer_subtree) {
		/* Slice the hyperrect to get the hyperrect of the nearer subtree */
		dummy = *nearer_hyperrect_coord;
		*nearer_hyperrect_coord = node->pos[dir];
		/* Recurse down into nearer subtree */
		kd_nearest_i(nearer_subtree, pos, result, result_dist_sq, rect);
		/* Undo the slice */
		*nearer_hyperrect_coord = dummy;
	}

	/* Check the distance of the point at the current node, compare it
	 * with our best so far */
	dist_sq = 0;
	for(i=0; i < rect->dim; i++) {
		dist_sq += SQ(node->pos[i] - pos[i]);
	}
	if (dist_sq < *result_dist_sq && dist_sq > 0) {
		*result = node;
		*result_dist_sq = dist_sq;
	}

	if (farther_subtree) {
		/* Get the hyperrect of the farther subtree */
		dummy = *farther_hyperrect_coord;
		*farther_hyperrect_coord = node->pos[dir];
		/* Check if we have to recurse down by calculating the closest
		 * point of the hyperrect and see if it's closer than our
		 * minimum distance in result_dist_sq. */
		if (hyperrect_dist_sq(rect, pos) < *result_dist_sq) {
			/* Recurse down into farther subtree */
			kd_nearest_i(farther_subtree, pos, result, result_dist_sq, rect);
		}
		/* Undo the slice on the hyperrect */
		*farther_hyperrect_coord = dummy;
	}
}

struct kdres *kd_nearest(struct kdtree *kd, const double *pos)
{
	struct kdhyperrect *rect;
	struct kdnode *result;
	struct kdres *rset;
	double dist_sq;
	int i;

	if (!kd) return 0;
	if (!kd->rect) return 0;

	/* Allocate result set */
	if(!(rset = malloc(sizeof *rset))) {
		return 0;
	}
	if(!(rset->rlist = alloc_resnode())) {
		free(rset);
		return 0;
	}
	rset->rlist->next = 0;
	rset->tree = kd;

	/* Duplicate the bounding hyperrectangle, we will work on the copy */
	if (!(rect = hyperrect_duplicate(kd->rect))) {
		kd_res_free(rset);
		return 0;
	}

	/* Our first guesstimate is the root node */
	result = kd->root;
	dist_sq = 0;
	for (i = 0; i < kd->dim; i++)
		dist_sq += SQ(result->pos[i] - pos[i]);

	/* Search for the nearest neighbour recursively */
	kd_nearest_i(kd->root, pos, &result, &dist_sq, rect);

	/* Free the copy of the hyperrect */
	hyperrect_free(rect);

	/* Store the result */
	if (result) {
		if (rlist_insert(rset->rlist, result, -1.0) == -1) {
			kd_res_free(rset);
			return 0;
		}
		rset->size = 1;
		kd_res_rewind(rset);
		return rset;
	} else {
		kd_res_free(rset);
		return 0;
	}
}


/* ---- nearest N search ---- */

struct kdres *kd_nearest_n(struct kdtree *kd, const double *pos, int num)
{   double range = 999999.0;
	int ret;
	struct kdres *rset;
    struct rheap *h;
    h = initHeap();
	if(!(rset = malloc(sizeof *rset))) {
		return 0;
	}
	if(!(rset->rlist = alloc_resnode())) {
		free(rset);
		return 0;
	}
	rset->rlist->next = 0;
	rset->tree = kd;
	if((ret = find_nearest_n(kd->root, pos, range, num, h, kd->dim)) == -1) {
		kd_res_free(rset);
		return 0;
	}
	rset->size = ret;
	kd_res_rewind(rset);
	return rset;
}

struct kdres *kd_nearest_range(struct kdtree *kd, const double *pos, double range)
{
	int ret;
	struct kdres *rset;

	if(!(rset = malloc(sizeof *rset))) {
		return 0;
	}
	if(!(rset->rlist = alloc_resnode())) {
		free(rset);
		return 0;
	}
	rset->rlist->next = 0;
	rset->tree = kd;

	if((ret = find_nearest(kd->root, pos, range, rset->rlist, 0, kd->dim)) == -1) {
		kd_res_free(rset);
		return 0;
	}
	rset->size = ret;
	kd_res_rewind(rset);
	return rset;
}

void kd_res_free(struct kdres *rset)
{
	clear_results(rset);
	free_resnode(rset->rlist);
	free(rset);
}

int kd_res_size(struct kdres *set)
{
	return (set->size);
}

void kd_res_rewind(struct kdres *rset)
{
	rset->riter = rset->rlist->next;
}

int kd_res_end(struct kdres *rset)
{
	return rset->riter == 0;
}

int kd_res_next(struct kdres *rset)
{
	rset->riter = rset->riter->next;
	return rset->riter != 0;
}

int *kd_res_item(struct kdres *rset, double *pos)
{
	if(rset->riter) {
		if(pos) {
			memcpy(pos, rset->riter->item->pos, rset->tree->dim * sizeof *pos);
		}
		return rset->riter->item->data;
	}
	return 0;
}

int *kd_res_itemf(struct kdres *rset, float *pos)
{
	if(rset->riter) {
		if(pos) {
			int i;
			for(i=0; i<rset->tree->dim; i++) {
				pos[i] = rset->riter->item->pos[i];
			}
		}
		return rset->riter->item->data;
	}
	return 0;
}

int *kd_res_item_data(struct kdres *set)
{
	return kd_res_item(set, 0);
}

/* ---- hyperrectangle helpers ---- */
static struct kdhyperrect* hyperrect_create(int dim, const double *min, const double *max)
{
	size_t size = dim * sizeof(double);
	struct kdhyperrect* rect = 0;

	if (!(rect = malloc(sizeof(struct kdhyperrect)))) {
		return 0;
	}

	rect->dim = dim;
	if (!(rect->min = malloc(size))) {
		free(rect);
		return 0;
	}
	if (!(rect->max = malloc(size))) {
		free(rect->min);
		free(rect);
		return 0;
	}
	memcpy(rect->min, min, size);
	memcpy(rect->max, max, size);

	return rect;
}

static void hyperrect_free(struct kdhyperrect *rect)
{
	free(rect->min);
	free(rect->max);
	free(rect);
}

static struct kdhyperrect* hyperrect_duplicate(const struct kdhyperrect *rect)
{
	return hyperrect_create(rect->dim, rect->min, rect->max);
}

static void hyperrect_extend(struct kdhyperrect *rect, const double *pos)
{
	int i;

	for (i=0; i < rect->dim; i++) {
		if (pos[i] < rect->min[i]) {
			rect->min[i] = pos[i];
		}
		if (pos[i] > rect->max[i]) {
			rect->max[i] = pos[i];
		}
	}
}

static double hyperrect_dist_sq(struct kdhyperrect *rect, const double *pos)
{
	int i;
	double result = 0;

	for (i=0; i < rect->dim; i++) {
		if (pos[i] < rect->min[i]) {
			result += SQ(rect->min[i] - pos[i]);
		} else if (pos[i] > rect->max[i]) {
			result += SQ(rect->max[i] - pos[i]);
		}
	}

	return result;
}

/* ---- static helpers ---- */

#ifdef USE_LIST_NODE_ALLOCATOR
/* special list node allocators. */
static struct res_node *free_nodes;

#ifndef NO_PTHREADS
static pthread_mutex_t alloc_mutex = PTHREAD_MUTEX_INITIALIZER;
#endif

static struct res_node *alloc_resnode(void)
{
	struct res_node *node;

#ifndef NO_PTHREADS
	pthread_mutex_lock(&alloc_mutex);
#endif

	if(!free_nodes) {
		node = malloc(sizeof *node);
	} else {
		node = free_nodes;
		free_nodes = free_nodes->next;
		node->next = 0;
	}

#ifndef NO_PTHREADS
	pthread_mutex_unlock(&alloc_mutex);
#endif

	return node;
}

static void free_resnode(struct res_node *node)
{
#ifndef NO_PTHREADS
	pthread_mutex_lock(&alloc_mutex);
#endif

	node->next = free_nodes;
	free_nodes = node;

#ifndef NO_PTHREADS
	pthread_mutex_unlock(&alloc_mutex);
#endif
}
#endif	/* list node allocator or not */


/* inserts the item. if dist_sq is >= 0, then do an ordered insert */
/* TODO make the ordering code use heapsort */
static int rlist_insert(struct res_node *list, struct kdnode *item, double dist_sq)
{
	struct res_node *rnode;

	if(!(rnode = alloc_resnode())) {
		return -1;
	}
	rnode->item = item;
	rnode->dist_sq = dist_sq;

	if(dist_sq >= 0.0) {
		while(list->next && list->next->dist_sq < dist_sq) {
			list = list->next;
		}
	}
	rnode->next = list->next;
	list->next = rnode;
	return 0;
}

static void clear_results(struct kdres *rset)
{
	struct res_node *tmp, *node = rset->rlist->next;

	while(node) {
		tmp = node;
		node = node->next;
		free_resnode(tmp);
	}

	rset->rlist->next = 0;
}




