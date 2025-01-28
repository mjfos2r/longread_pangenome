class CursedDict(defaultdict):
    """
    may the divine forgive me for such an egregious hack.
    Basically to keep from refactoring everything via a custom factory we just
    commit this atrocity and instantiate two hard coded keys and then treat the
    rest as a defaultdict(list). 
    This is going to break many thing probably but I hope it fixes at least two problems.
    """
    def __init__(self):
        super().__init__(list) # super().__scuffed__
        self["assembly_id"] = ""
        self["contig_length"] = int()

    def __getitem__(self, key):
        #if key in ("assembly_id", "contig_length"):
        #    return super().__getitem__(key) # no hardcoding required?
        return super().__getitem__(key)

    def __setitem__(self, key, value):
        #this is so bad lmao
        if key == "assembly_id" and not isinstance(value, str):
            raise ValueError("assembly_id must be a string")
        elif key == "contig_length" and not isinstance(value, int):
            raise ValueError("contig_length must be an integer")
        super().__setitem__(key, value)
