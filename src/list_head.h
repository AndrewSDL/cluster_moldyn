#ifndef _list_head_H
#define	_list_head_H

struct list_head {
    struct list_head *next, *prev;
    void* data;
};


#endif	/* _list_head_H */

