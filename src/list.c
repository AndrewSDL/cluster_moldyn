#include "list.h"

// - implementeaza operatii elementare cu liste santinela circulare dubluinlantuite,
// in care informatia utila este indicata de membrul 'data'

struct list_head* new_list_head()
{
	struct list_head* list = (struct list_head*) malloc( sizeof( struct list_head ) );
	init_list_head( list );
	return list;
}

struct list_head* new_list_elem(void* data)
{
	struct list_head* elem = (struct list_head*) malloc( sizeof( struct list_head ) );
	elem->next = elem;
        elem->prev = elem;
        elem->data = data;
	return elem;
}

void list_add_between(struct list_head *new_lh,
                              struct list_head *prev,
                              struct list_head *next)
{
        next->prev = new_lh;
        new_lh->next = next;
        new_lh->prev = prev;
        prev->next = new_lh;
}

void list_add_first(struct list_head *new_lh, struct list_head *head)
{
        list_add_between(new_lh, head, head->next);
}

void list_add_last(struct list_head *new_lh, struct list_head *head)
{
        list_add_between(new_lh, head->prev, head);
}

void list_add_data_first(void* data, struct list_head *head)
{
	struct list_head* elem = new_list_elem( data );
	list_add_first(elem, head);
}



// copie in new_list pe old_list
void list_copy(struct list_head* new_list, struct list_head*
old_list )
{
	struct list_head *i;

	list_head_for_each( i, old_list )
	{
		list_add_data_last( i->data, new_list );
	}
}


// copie in new_list pe old_list folosind functia de copiere copy_data pentru copierea campului data
void list_copy2(struct list_head* new_list, struct list_head* old_list, void* (*copy_data)(void*) )
{
	struct list_head *i;

	list_head_for_each( i, old_list )
	{
		list_add_data_last( copy_data(i->data), new_list );
	}
}

//  intoarce elementul listei pentru care functia de comparare returneaza egalitate,
struct list_head* list_find( void* data, struct list_head* head, int (*cmp)(void*,void*) )
{
	struct list_head* i;

	list_head_for_each( i, head )
	{
		if( cmp( data, i->data ) == 0 )
			return i;
	}

	return NULL;
}

//  intoarce pointerul catre membrul data al elementului listei pentru care
// functia de comparare cu 'data' returneaza egalitate,
void* list_find_data( void* data, struct list_head* head, int (*cmp)(void*,void*) )
{
	//printf( "list_find_data \n" );

	struct list_head* elem = list_find( data, head, cmp );

	//printf( "list_find_data end \n" );

	if( elem )
		return elem->data;
	return NULL;
}


void __list_del(struct list_head * prev, struct list_head * next)
{
        next->prev = prev;
        prev->next = next;
}

void list_del2( struct list_head* entry )
{
    list_del(entry);
    free(entry);
}

void list_destroy( struct list_head* head )
{
	struct list_head *i, *n;
	list_head_for_each_safe( i, n, head )
	{
		list_del(i);
		free( i );
	}
}

void free_list( struct list_head* head )
{
	list_destroy( head );
	free( head );
}



void list_print( struct list_head* head, char* (*to_string)(void*), char* sep )
{
	struct list_head* i;
	list_head_for_each( i, head )
		printf("%s%s", to_string(i->data), sep);
	printf("\n");
}

int list_length( struct list_head* head )
{
	int k = 0;
	struct list_head *i;

	list_head_for_each( i, head )
		k++;
	return k;
}
