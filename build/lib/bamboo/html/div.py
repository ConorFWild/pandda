

class Div(object):


    def __init__(self, classes):
        self.classes = classes
        self.children = []

    def add_class(self, class_name):
        self.classes.append(class_name)

    def add(self, other):
        """Create a div object and add it to the children"""
        self.children.append(other)

    def format(self):
        return "<div class='{}'>\n".format(' '.join(self.classes))+\
                '  '+'\n'.join([c.format() for c in self.children]).replace('\n','\n  ')+\
                "\n</div>"


class Col(Div):

    def set_width(self, size='lg', value=12):
        self.classes.append('col-{}-{}'.format(size.strip(), str(value).strip()))

class Row(Div):
    pass

