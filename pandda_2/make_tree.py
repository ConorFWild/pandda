class MakeTree:
    def __init__(self, overwrite=False):
        self.overwrite = overwrite

    def __call__(self, tree):
        tree.make(overwrite=self.overwrite)