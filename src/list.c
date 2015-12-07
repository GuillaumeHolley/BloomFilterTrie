#include "./../lib/list.h"

List *List_create()
{
    return calloc(1, sizeof(List));
}

void List_destroy(List *list)
{
    LIST_FOREACH(list, first, next, cur) {
        if(cur->prev) {
            free(cur->prev);
        }
    }

    free(list->last);
    free(list);
}


void List_clear(List *list)
{
    LIST_FOREACH(list, first, next, cur) {
        free(cur->value);
    }
}


void List_clear_destroy(List *list)
{
    List_clear(list);
    List_destroy(list);
}


void List_push(List *list, void *value)
{
    ListNode *node = calloc(1, sizeof(ListNode));
    ASSERT_NULL_PTR(node, "List_push()\n")

    node->value = value;

    if(list->last == NULL) {
        list->first = node;
        list->last = node;
    } else {
        list->last->next = node;
        node->prev = list->last;
        list->last = node;
    }

    list->count++;

    return;
}

void *List_remove_first(List *list)
{
    ListNode* node = list->first;
    return node != NULL ? List_remove(list, node) : NULL;
}

void *List_pop(List *list)
{
    ListNode *node = list->last;
    return node != NULL ? List_remove(list, node) : NULL;
}

void *List_remove(List *list, ListNode *node)
{
    void *result = NULL;

    if ((list->first == NULL) && (list->last == NULL)) ERROR("List is empty.")
    if (node == NULL) ERROR("node can't be NULL")

    if(node == list->first && node == list->last) {
        list->first = NULL;
        list->last = NULL;
    } else if(node == list->first) {
        list->first = node->next;
        if (list->first == NULL) ERROR("Invalid list, somehow got a first that is NULL.")
        list->first->prev = NULL;
    } else if (node == list->last) {
        list->last = node->prev;
        if (list->last == NULL) ERROR("Invalid list, somehow got a next that is NULL.")
        list->last->next = NULL;
    } else {
        ListNode *after = node->next;
        ListNode *before = node->prev;
        after->prev = before;
        before->next = after;
    }

    list->count--;
    result = node->value;
    free(node);

    return result;
}

void *List_first_become_last(List *list)
{
    ListNode* node = list->first;

    if ((list->first == NULL) && (list->last == NULL)) ERROR("List is empty.")
    if (node == NULL) ERROR("node can't be NULL")

    if (node != list->last){

        list->first = node->next;
        if (list->first == NULL) ERROR("Invalid list, somehow got a first that is NULL.")
        list->first->prev = NULL;

        list->last->next = node;
        node->prev = list->last;
        node->next = NULL;
        list->last = node;
    }

    return node->value;
}
