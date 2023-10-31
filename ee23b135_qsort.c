/**
 * EE23B135 Kaushik G Iyer
 * 05/09/2023
 * 
 * Sorts a doubly linked list (which is read from a file) using the quicksort algorithm
 * 
 * Expected input:
 *  Just one argument containing the path to the file with integers
 *  NOTE: The file is expected to have whitespace seperated integers
 *
 * Outputs:
 *  Space separated list of sorted integers (Ascending)
 * 
*/ 
#include <stdio.h>
#include <stdlib.h>

struct Node{
    int val;
    struct Node* prev;
    struct Node* next;
};
struct Node* NODE__init(int val);
void NODE__destroy(struct Node* head);
struct Node* NODE__adopt(struct Node* parent, struct Node* child);
struct Node* NODE__disownChild(struct Node* parent);
void NODE__insertRight(struct Node* parent, struct Node* node);
void NODE__insertLeft(struct Node* parent, struct Node* node);
void NODE__print(struct Node* head);
struct Node* NODE__last(struct Node* head);


/**
 * Basically since when we call partition, the order changes, we must return the new head, new end also
*/
struct PartitionNodeTuple{ // Just a struct that contains 3 pointers to nodes
    struct Node* newHead;
    struct Node* newEnd;
    struct Node* partition;
};
struct PartitionNodeTuple partition(struct Node* start, struct Node* end);
struct Node* quicksort(struct Node* start, struct Node* end);

int main(int argc, char** argv){
    if (argc < 2){
        printf("Please give file name vro\n");
        return 1;
    }
    FILE* file = fopen(argv[1], "r");
    if (file == NULL){
        printf("Please give valid file name vro\n");
        return 1;
    }

    struct Node* head = NODE__init(-1); // Just for our convenience we have a dummy node as our head

    int num;
    struct Node* tail = head;
    while (fscanf(file, "%d ", &num) == 1){
        tail = NODE__adopt(tail, NODE__init(num));
    }
    fclose(file);   

    if (head->next == NULL){ // i.e. We could not read any integers
        printf("Please provide a file with whitespace separated integers only!\n");
        return 1;
    }

    // NOTE: The quicksort function returns the pointer to the new start (so if you wanted to remove the dummy node completely, you would have to say head = quicksort(head, last) as a new node may become the first node)
    quicksort(head->next, tail); // The actual linked list that we took as input starts from head.next
    NODE__print(head->next);

    NODE__destroy(head); // Just free up all the memory we allocated
    return 0;
};


/**
 * Basically just manipulates the array such that all numbers smaller than the pivot are to the left of it and all number bigger are to the right
 * @returns Reference to the pivot 
 * 
 * NOTE: It is expected that end is a desendant of start (i.e. start(->next) finite amount of times == end)
*/
struct PartitionNodeTuple partition(struct Node* start, struct Node* end){
    if (start == end){ // Handle the edge case where our list is only of 1 node
        struct PartitionNodeTuple pnt = {start, start, start};
        return pnt;
    }

    struct Node* pivot = start; // Our pivot node is the original first node (because it is convenient)
    struct Node* curr = start->next; // Keeps track of which node we are checking
    struct Node* tmp;
    while (curr != end){ // We check all nodes except for the last node, that one is special
        tmp = curr->next;
        if (curr->val < pivot->val){ // This means we have to move `curr` to the left of `pivot`
            
            // Ok so basically let us just keep inserting smaller nodes to the left of start...
            // That should always ensure that the newly inserted node is the new start
            NODE__disownChild(curr->prev); // Pop out curr basically
            NODE__insertLeft(start, curr);
            start = curr; // Since we added it to the left of start, it is now the new start
        }
        // Note that if that value of the node we are seeing is greater than the pivot value we don't need to to anything
        // (Since it is already to the right of the pivot)
        curr = tmp;
    }

    // The end is special since if we move it we must update the new end value
    if (end->val < pivot->val){
        tmp = end->prev; // Keep track of the new end
        NODE__disownChild(end->prev); // Pop out end
        NODE__insertLeft(start, end); // Add it to the left of pivot
        start = end;
        end = tmp;
    }

    struct PartitionNodeTuple pnt = {start, end, pivot};
    return pnt;
}


/**
 * Sorts the list using quicksort bruv
 * NOTE: It is expected that end is a desendant of start (i.e. start(->next) finite amount of times == end)
 * 
 * @returns Reference to the new head of the list
*/
struct Node* quicksort(struct Node* start, struct Node* end){
    if (start == end){return start;} // A one element list is already sorted :)
    
    struct PartitionNodeTuple pnt = partition(start, end);
    
    if (pnt.newEnd != pnt.partition)
        quicksort(pnt.partition->next, pnt.newEnd); // Sort the right of the list
    // If our partition is at the end, there are no numbers to the right of it (i.e. nothing to be sorted)
    
    
    if (pnt.newHead != pnt.partition)
        return quicksort(pnt.newHead, pnt.partition->prev); // Return the new head after sorting the left of the list

    // If our partition is at the beginning, there are no numbers to the left of it (i.e. nothing to be sorted)
    return pnt.newHead;
}


/**
 * Allocate memory and initialize the node
 * Q. Why not just return a copy of the new node (return struct Node) ?
 * A. For some reason this led to previously created nodes being overwritten. How? idk, I was probably doing something stupid.
*/
struct Node* NODE__init(int val){
    struct Node* node = malloc(sizeof(struct Node));
    node->next = NULL;
    node->prev = NULL;
    node->val = val;
    return node;
}


/**
 * Frees up all the allocated memory (Since we are creating nodes using malloc)
*/
void NODE__destroy(struct Node* head){
    struct Node* next;
    while (head != NULL){
        next = head->next;
        free(head);
        head = next;
    }    
}


/**
 * Sets the `child` as the next of `parent`
 * And sets the `parent` as prev of `child`
*/
struct Node* NODE__adopt(struct Node* parent, struct Node* child){
    if (parent != NULL) parent->next = child;
    if (child != NULL) child->prev = parent;
    return child;
}


/**
 * Basically pops out the next node in the list
*/
struct Node* NODE__disownChild(struct Node* parent){
    struct Node* child = parent->next;
    if (child == NULL) return child; // i.e there is no child to disown :(

    NODE__adopt(parent, child->next);
    
    child->prev = NULL;
    child->next = NULL;

    return child;
}


/**
 * Inserts the insertee just after the node
*/
void NODE__insertRight(struct Node* node, struct Node* insertee){
    struct Node* next = node->next;
    NODE__adopt(node, insertee);
    NODE__adopt(insertee, next);
}


/**
 * Inserts the insertee just before the node
*/
void NODE__insertLeft(struct Node* node, struct Node* insertee){
    struct Node* prev = node->prev;
    NODE__adopt(insertee, node);
    NODE__adopt(prev, insertee);
}


/**
 * Prints all the values from the list (Space separated)
*/
void NODE__print(struct Node* head){
    if (head == NULL){
        printf("\n");
        return;
    }
    printf("%d ", head->val);
    NODE__print(head->next);
}


/**
 * Finds the last node in the list
*/
struct Node* NODE__last(struct Node* head){
    struct Node* curr = head;
    while (curr->next != NULL) curr = curr->next;
    return curr;
}