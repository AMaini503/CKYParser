class Tree:

    """Binary tree implementation
    
    The leaves are the nodes with left and right == None

    """

    def __init__(self, val):
        self.val = val
        self.left = None
        self.right = None
        
        # `down` is treated separately to distinguish between internal nodes (non-termianls)
        # and leaves (words)
        # This is helpful when converting this tree into array of nodes
        self.down = None
    def SetLeftChild(self, left):
        assert(isinstance(left, Tree))
        self.left = left

    def SetRightChild(self, right):
        assert(isinstance(right, Tree))
        self.right = right       

    def SetDownChild(self, down):
        # The leaves (words) are always strings
        assert(isinstance(down, str)) 

        # There should be no left and right child
        assert(self.left is None)
        assert(self.right is None)

        self.down = down

    def isFringe(self):
        if(self.down == None):
            return False
        else:
            return True

    def getVal(self):
        return self.val
    
    def toArray(self):
        """Converts tree structure to an array"""

        # If this node is a Fringe, its child is a word
        if(self.isFringe() == True):
            return [self.val, self.down]
        # Every other node has 2 children
        else:
            return [self.val, self.left.toArray(), self.right.toArray()]
