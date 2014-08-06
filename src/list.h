#ifndef _list_H
#define	_list_H

#ifdef	__cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>

#include "list_head.h"

// - implementeaza operatii elementare cu liste santinela circulare dubluinlantuite,
// in care informatia utila este indicata de membrul 'data'


#define LIST_HEAD_INIT(name) { &(name), &(name) , NULL }

#define init_list_head(ptr) do { \
        (ptr)->next = (ptr); (ptr)->prev = (ptr); (ptr)->data = NULL;\
} while (0)

struct list_head* new_list_head();
struct list_head* new_list_elem(void* data);

void list_add_between(struct list_head *new_lh, struct list_head *prev, struct list_head *next);
void list_add_first(struct list_head *new_lh, struct list_head *head);
void list_add_last(struct list_head *new_lh, struct list_head *head);

void list_add_data_first(void* data, struct list_head *head);

// adauga un element la sf listei, cu campul 'data'
#define list_add_data_last( lst_data, lst_head) \
{   \
	struct list_head* elem = (struct list_head*) malloc( sizeof( struct list_head ) );  \
        struct list_head* _prev = (lst_head)->prev;   \
	elem->data = lst_data;  \
	(lst_head)->prev = elem;  \
        elem->next = (lst_head);  \
        elem->prev = _prev;     \
        _prev->next = elem;     \
}

void list_copy(struct list_head* new_list, struct list_head* old_list );
void list_copy2(struct list_head* new_list, struct list_head* old_list, void* (*copy_data)(void*) );

struct list_head* list_find( void* data, struct list_head* head, int (*cmp)(void*,void*) );
void* list_find_data( void* data, struct list_head* head, int (*cmp)(void*,void*) );

void __list_del(struct list_head * prev, struct list_head * next);

#define list_del(entry)        { (entry)->next->prev = (entry)->prev; (entry)->prev->next = (entry)->next; }
void list_del2( struct list_head* entry );

#define list_empty(head)  ((head)->next == (head))

void list_destroy( struct list_head* head );
void free_list( struct list_head* head );

void list_print( struct list_head* head, char* (*to_string)(void*),char* sep );
int list_length( struct list_head* head );


#define list_head_for_each(pos, head) \
        for (pos = (head)->next; pos != (head); pos = pos->next)

// parcurge lista in sens invers
#define list_head_for_each_prev(pos, head) \
        for (pos = (head)->prev; pos != (head); pos = pos->prev)

#define list_head_for_each_safe(pos, n, head) \
        for (pos = (head)->next, n = pos->next; pos != (head); \
                pos = n, n = pos->next)

// parcurge 2 liste in paralel
#define list_for_each2(pos1, pos2, head1, head2) \
        for (pos1 = (head1)->next, pos2 = (head2)->next; pos1 != (head1) && pos2 != (head2); pos1 = pos1->next, pos2 = pos2->next)

// parcurge cele 2 liste in paralel de la coada la cap
#define list_for_each2_prev(pos1, pos2, head1, head2) \
        for (pos1 = (head1)->prev, pos2 = (head2)->prev; pos1 != (head1) && pos2 != (head2); pos1 = pos1->prev, pos2 = pos2->prev)

#define list_for_each_safe2(pos1, pos2, n1, n2, head1, head2) \
        for (pos1 = (head1)->next, n1 = pos1->next, pos2 = (head2)->next, n2 = pos2->next; pos1 != (head1) && pos2 != (head2); \
               pos1 = n1, n1 = pos1->next, pos2 = n2, n2 = pos2->next)





#ifdef	__cplusplus
}
#endif

#endif	/* _list_H */

