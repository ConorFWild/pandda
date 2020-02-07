import os, sys
import base64, hashlib
import socket
import threading

# CONNECTION CLASSES

class ConnectionReferenceObject(object):
    """ Object for storing socket information """
    def __init__(self):
        self.sock = None
        self.conn = None
        self.lock = None
        self.addr = None

class SocketReferenceObject(object):
    """ Object for passing Incoming and Outgoing Socket Objects """
    def __init__(self):
        """Spaces for incoming/outgoing sockets, connections, locks"""
        self.incoming = ConnectionReferenceObject()
        self.outgoing = ConnectionReferenceObject()

class ThreadWithReturnValue(threading.Thread):
    """ Define Thread so that it has a return value """
    def __init__(self, group=None, target=None, name=None,
                 args=(), kwargs={}, verbose=False):
        threading.Thread.__init__(self, group, target, name, args, kwargs, verbose)
        self._return = None
    def run(self):
        if self._Thread__target is not None:
            self._return = self._Thread__target(*self._Thread__args,
                                                **self._Thread__kwargs)
    def join(self):
        threading.Thread.join(self)
        return self._return

# CONNECTION FUNCTIONS

def create_handshake_resp(handshake):
    """ Respond to handshake to establish connection """
    final_line = ""
    lines = handshake.splitlines()
    for line in lines:
        parts = line.partition(": ")
        if parts[0] == "Sec-WebSocket-Key":
            key = parts[2]

    magic = '258EAFA5-E914-47DA-95CA-C5AB0DC85B11'
    accept_key = base64.b64encode(hashlib.sha1(key+magic).digest())

    return (
        "HTTP/1.1 101 Switching Protocols\r\n"
        "Upgrade: WebSocket\r\n"
        "Connection: Upgrade\r\n"
        "Sec-WebSocket-Accept: " + accept_key + "\r\n\r\n")

# DECODING FUNCTIONS

def decode_message(string_stream_in):
    """ Decode a WebSocket data frame into a string """

    #turn string values into opererable numeric byte values
    byteArray = [ord(character) for character in string_stream_in]
    datalength = byteArray[1] & 127
    indexFirstMask = 2
    if datalength == 126:
        indexFirstMask = 4
    elif datalength == 127:
        indexFirstMask = 10
    masks = [m for m in byteArray[indexFirstMask : indexFirstMask+4]]
    indexFirstDataByte = indexFirstMask + 4
    decoded_chars = []
    i = indexFirstDataByte
    j = 0
    while i < len(byteArray):
        decoded_chars.append( chr(byteArray[i] ^ masks[j % 4]) )
        i += 1
        j += 1

    return ''.join(decoded_chars)

def encode_message_to_send(raw_data):
    """Turns a string into a WebSocket data frame. Returns a bytes(). 'data' is a string"""

    #determine the size of the data we were told to send
    data_length = len(raw_data)
    output_bytes = bytearray()
    output_bytes.append(0x81) #0x81 = text data type
    if data_length < 0x7D:
        #a nice short length
        output_bytes.append(len(raw_data))
    elif data_length >= 0x7E and len(raw_data) < 0xFFFF:
        #two additional bytes of length needed
        output_bytes.append(0x7E)
        output_bytes.append(data_length >> 8 & 0xFF)
        output_bytes.append(data_length & 0xFF)
    else:
        #eight additional bytes of length needed
        output_bytes.append(0x7F)
        output_bytes.append(data_length >> 56 & 0xFF)
        output_bytes.append(data_length >> 48 & 0xFF)
        output_bytes.append(data_length >> 40 & 0xFF)
        output_bytes.append(data_length >> 32 & 0xFF)
        output_bytes.append(data_length >> 24 & 0xFF)
        output_bytes.append(data_length >> 16 & 0xFF)
        output_bytes.append(data_length >> 8 & 0xFF)
        output_bytes.append(data_length & 0xFF)
    #tack on the raw data now
    for byte in raw_data:
        output_bytes.append(ord(byte))
    return bytes(output_bytes)

def send_message_to_connection(message, connection):
    """ Message is a String. Connection is a ConnectionReferenceObject. """
    msg = encode_message_to_send(message)
    connection.lock.acquire()
    connection.conn.send(msg)
    connection.lock.release()

def establish_outgoing_connection(outgoing_connection, ADDR, PORT):
    """ Establish a connection with COOT """

    print "Outgoing Connection: ", ADDR, ' on ', PORT
    print "WAITING FOR OUTGOING CONNECTION TO BE MADE"

    outgoing_connection.sock = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
    outgoing_connection.sock.setsockopt(socket.SOL_SOCKET,socket.SO_REUSEADDR,1)
    outgoing_connection.sock.bind((ADDR,PORT))
    outgoing_connection.sock.listen(1)

    # Handle Request
    outgoing_connection.conn, outgoing_connection.addr = outgoing_connection.sock.accept()
    print 'Outgoing Connection made to: ', outgoing_connection.addr
    data = outgoing_connection.conn.recv(1024)
    print 'Received Msg: ', data

    return outgoing_connection

def establish_incoming_connection(incoming_connection, ADDR, PORT):
    """ Establish a connection with the web client """

    print "Incoming Connection: ", ADDR, ' on ', PORT
    print "WAITING FOR INCOMING CONNECTION TO BE MADE"

    incoming_connection.sock = socket.socket()
    incoming_connection.sock.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
    incoming_connection.sock.bind((ADDR, PORT))
    incoming_connection.sock.listen(1)

    # Handle Request
    incoming_connection.conn, incoming_connection.addr = incoming_connection.sock.accept()
    print 'Incoming Connection made to: ', incoming_connection.addr
    data = incoming_connection.conn.recv(1024)
    print 'Received Msg: ', data
    response = create_handshake_resp(data)
    incoming_connection.conn.sendto(response, incoming_connection.addr)

    # Create a lock object
    incoming_connection.lock = threading.Lock()

    return incoming_connection

def close_outgoing_connection(outgoing_connection):
    """ Close Outgoing Connection """

    outgoing_connection.conn.shutdown(2)
    outgoing_connection.sock.close()

    return 0

def close_incoming_connection(incoming_connection):
    """ Close Incoming Connection """

    print 'Closing Client:', incoming_connection.addr

    incoming_connection.lock.acquire()
    incoming_connection.addr = ''
    incoming_connection.sock.close()
    incoming_connection.lock.release()

    return 0








