class node:
    def __init__(self):
        self.data = None # contains the data
        self.next = None # contains the reference to the next node

class linked_list:
    def __init__(self):
        self.head = None

    def insert_node(self, data):
        new_node = node() # create a new node
        new_node.data = data
        new_node.next = self.head # link the new node to the 'previous' node.
        self.head = new_node #  set the current node to the new one.

    def list_print(self):
        node = self.head
        while node:
            print node.data
            node = node.next
