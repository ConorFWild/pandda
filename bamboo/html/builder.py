
import jinja2

import bamboo.html.section

class HtmlPage(objects):
    """Class for programatically generating an HTML page"""

    def __init__(self, title=None):

        self.header = bamboo.html.section.Header(title=title)
        self.content = bamboo.html.section.Content()
        self.footer = bamboo.html.section.Footer()

    def format(self):

        return '<html>\n'+\
                self.header.format()+\
                self.content.format()+\
                self.footer.format()+\
               '</html>'
