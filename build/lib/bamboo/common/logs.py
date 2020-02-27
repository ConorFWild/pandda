import os, sys, datetime


class Bar(object):


    def __init__(self, width=40, body='-', head='>>>'):
        self.width = width
        self.body  = body
        self.head  = head

    def __call__(self, width=None, body=None, head=None):
        # Allow defaluts to be overridden
        if width is None: width = self.width
        if body is None:  body  = self.body
        if head is None:  head  = self.head
        # Modify width to allow for head
        width = width - len(head)
        if width < 1: width = 0
        return body*width + head


class Heading(object):


    def __init__(self, width=100, spacer='#', decorator=' <~~~> ', side_width=3):
        """Spacer for printing neatly formatted headings"""

        assert width > 10
        self._decorator     = decorator
        self._spacer        = spacer
        self._width         = width
        self._side_width    = side_width
        self._content_width = width-(2*self._side_width)

        self._side_padding  = self._side_width*spacer

    def __call__(self, text='', spacer=False, blank=False):

        # Allow for need to override defaults
        content_width = max(len(text)+2, self._content_width)
        total_width = content_width + (2*self._side_width)

        out = [self._text_line(text=text, text_width=content_width)]
        if spacer:
            out.insert(0, self._inner_line(text_width=content_width))
            out.append(out[0])
        out.insert(0, self._outer_line(width=total_width))
        out.append(out[0])
        if blank:
            out.insert(0, self._blank_line())
            out.append(out[0])
        return '\n'.join(out)

    def _blank_line(self):
        return ''

    def _text_line(self, text, text_width):
        return self._side_padding + text.center(text_width, ' ') + self._side_padding

    def _inner_line(self, text_width):
        return self._text_line(text='', text_width=text_width)

    def _outer_line(self, width):
        return self._decorator.center(width, self._spacer)


class Log(object):


    def __init__(self, log_file=None, stdout=sys.stdout, verbose=False, silent=False):
        """Log Object for writing to logs...?"""
        assert not (silent and verbose), 'cannot be both silent and verbose...'
        # File that the log will be written to
        self.log_file = log_file
        self.verbose = verbose
        self.silent = silent

        # Set default heading
        self.set_bar_and_heading(bar        = Bar(),
                                 heading    = Heading(spacer='#', decorator=' <~~~> '),
                                 subheading = Heading(spacer='-', decorator=' *** ')
                                )

    def __call__(self, message, show=True):
        """Log message to file, and mirror to stdout if not silent and verbose or show is True"""
        message = str(message)+'\n'
        if (not self.silent) and (self.verbose or show):
            self.show(message)
        self.write(message=message, mode='a')

    def delete(self):
        """Delete the current log file"""
        if self.log_file and os.path.exists(self.log_file):
            os.remove(self.log_file)
        return self

    def heading(self, message, spacer=False, blank=True):
        """Print message as a heading/divider"""
        text = self._heading(text=str(message), spacer=spacer, blank=blank)
        self(text)

    def subheading(self, message, spacer=False, blank=True):
        """Print message as a heading/divider"""
        text = self._subheading(text=str(message), spacer=spacer, blank=blank)
        self(text)

    def bar(self, blank_before=False, blank_after=False, width=None, body=None, head=None):
        """Print divider bar"""
        text = '\n'*blank_before + self._bar(width=width, body=body, head=head) + '\n'*blank_after
        self(text)

    def set_bar_and_heading(self, bar=None, heading=None, subheading=None):
        if bar is not None:
            assert isinstance(bar, Bar)
            self._bar = bar
        if heading is not None:
            assert isinstance(heading, Heading)
            self._heading = heading
        if subheading is not None:
            assert isinstance(subheading, Heading)
            self._subheading = subheading
        return self

    def show(self, message):
        sys.stdout.write(message)

    def write(self, message, mode='a'):
        if not self.log_file: return
        message = message.replace('\r','')
        with open(self.log_file, mode) as fh: fh.write(message)

    def read_all(self):
        return open(self.log_file, 'r').read()

