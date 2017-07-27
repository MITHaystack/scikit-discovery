'''
This is a modified version of rforward.py, which can be found at:
https://github.com/paramiko/paramiko/blob/25dd096da065b1bc2f35c1a62d8a7055b022818b/demos/rforward.py

This version removes the command line interface, global variables, and creates tunnels in a non blocking thread.
An additional class has been added to facilitate the creation of tunnels

The original version is Copyright (C) 2008 Robey Pointer <robeypointer@gmail.com>

Since the original version is LGPL v2.1 or later, this is also LGPL v2.1 or later

Changes made to original

2017-03-30: Converted example into library for creating reverse ssh tunnels

2017-04-11: Updates to verbosity option and reverse tunneling

2017-07-10: Updates to documentation
'''


import socket
import threading
import select
import paramiko

def print_verbose(s, verbose=False):
    '''
    Print statement if verbose is True

    @param s: Statement to print
    @param verbose: Print only if verbose is True
    '''
    if verbose:
        print(s)


def handler(chan, host, port, verbose=False):
    '''
    Handler is responsible for sending and receiving data through ssh tunnel

    @param chan: SSH Channel for transferring data
    @param host: Address of remote host 
    @param port: Port to forward
    @param verbose: Print status information
    '''
    sock = socket.socket()
    try:
        sock.connect((host, port))
    except Exception as e:
        print_verbose('Forwarding request to %s:%d failed: %r' % (host, port, e), verbose)
        return
    
    print_verbose('Connected!  Tunnel open %r -> %r -> %r' % (chan.origin_addr,
                                                        chan.getpeername(), (host, port)), verbose)
    while True:
        r, w, x = select.select([sock, chan], [], [])
        if sock in r:
            data = sock.recv(1024)
            if len(data) == 0:
                break
            chan.send(data)
        if chan in r:
            data = chan.recv(1024)
            if len(data) == 0:
                break
            sock.send(data)
            
    chan.close()
    sock.close()
    print_verbose('Tunnel closed from %r' % (chan.origin_addr,), verbose)


def reverse_forward_tunnel(server_port, remote_host, remote_port, transport, check=30, verbose=False):
    '''
    Creates a reverse ssh tunnel

    @param server_port: Port on local host
    @param remote_host: Address of remote host
    @param remote_port: Port of remote host
    @param transport: SSH Transport
    @param check: Amount of time to wait in seconds when opening up a channel
    @param verbose: Print status information

    @return Thread running reverse ssh tunnel, event used to close ssh tunnel, 
            list of child threads started by main thread
    '''

    transport.request_port_forward('', server_port)
    
    event = threading.Event()
    child_threads = []
    
    def accept_tunnels(event):
        '''
        This function spawns new connections as they are needed
        
        @param event: When this event is set, this function will complete
        '''
        while not event.is_set():
            chan = transport.accept(check)
            if chan is None:
                continue
            thr = threading.Thread(target=handler, args=(chan, remote_host, remote_port, verbose))
            thr.setDaemon(True)
            thr.start()
            child_threads.append(thr)
                
    accept_thread = threading.Thread(target=accept_tunnels,args=[event])
    accept_thread.setDaemon(True)
    accept_thread.start()
    
    return accept_thread, event, child_threads
            

class ReverseTunnel(object):
    '''
    Create a reverse ssh tunnel
    '''

    def __init__(self, server_address, username, key_filename, server_port,
                 remote_host, remote_port, check=30, verbose=False):
        '''
        Initialize ReverseTunnel object

        @param server_address: Local server address
        @param username: Valid username on remote host
        @param key_filename: Filename of ssh key associated with remote host
        @param server_port: Local port
        @param remote_host: Address of remote host
        @param remote_port: Remote port
        @param check: Amount of time to wait in seconds when opening up a channel
        @param verbose: Print status information
        '''

        self.server_address = server_address
        self.username       = username
        self.key_filename   = key_filename
        self.server_port    = server_port
        self.remote_host    = remote_host
        self.remote_port    = remote_port
        self.check          = check
        self.verbose        = verbose

        self.ssh = None
        self.event = None

    
    def create_reverse_tunnel(self):
        '''
        Create the reverse tunnel
        '''
        self.ssh = paramiko.SSHClient()
        self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        self.ssh.connect(self.server_address, username=self.username, key_filename=self.key_filename)
        self.running_thread, self.event, self.child_threads = reverse_forward_tunnel(self.server_port, self.remote_host,
                                                                                     self.remote_port, self.ssh.get_transport(),
                                                                                     self.check)


    def __del__(self):
        '''
        Deconstructor
        '''
        if self.ssh != None:
            self.ssh.close()

        if self.event != None:
            self.event.set()
