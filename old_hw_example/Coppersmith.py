#!/usr/bin/sage

from sage.all import *
import struct
import re
from Crypto.PublicKey import RSA
from Crypto.Cipher import AES
from Crypto import Random

def parse_mpi(s,index):
    length = struct.unpack('I',s[index:index+4])[0]
    z = Integer(s[index+4:index+4+length].encode('hex'),16)
    return z, index+4+length


# Our "MPI" format consists of 4-byte integer length l followed by l bytes of binary key
def int_to_mpi(z):
    s = int_to_binary(z)
    return struct.pack('I',len(s))+s

# Horrible hack to get binary representation of arbitrary-length long int
def int_to_binary(z):
    s = ("%x"%z); s = (('0'*(len(s)%2))+s).decode('hex')
    return s

encrypt_header = '-----BEGIN PRETTY BAD ENCRYPTED MESSAGE-----\n'
encrypt_footer = '-----END PRETTY BAD ENCRYPTED MESSAGE-----\n'

# PKCS 7 pad message.
def pad(s,blocksize=AES.block_size):
    n = blocksize-(len(s)%blocksize)
    return s+chr(n)*n

# Encrypt string s using RSA encryption with AES in CBC mode.
# Generate a 256-bit symmetric key, encrypt it using RSA with PKCS1v1.5 type 1 padding, and prepend the MPI-encoded RSA ciphertext to the AES-encrypted ciphertext of the message.
def encrypt(rsakey,s):
    m = ZZ.random_element(2**256)

    k = ceil(rsakey.size()/8)-3-32
    EB = '0001' + 'ff'*k + '00' + "%x"%m

    output = int_to_mpi(rsakey.encrypt(int(EB,16),None))

    aeskey = int_to_binary(m)
    iv = Random.new().read(AES.block_size)
    cipher = AES.new(aeskey, AES.MODE_CBC, iv)

    output += iv + cipher.encrypt(pad(s))

    return encrypt_header + output.encode('base64') + encrypt_footer

"""
  I solved this problem by using Coppersmith's method. The key m takes up only
  256 bits of a key of size 2048 bits, and it is less than N^(1/3). Thus coppersmith's
  method should ideally work. Additionally, as Nadia's math shows, her lattice of
  dimension 4 works for essentially any m which is less has less than 1/6 of the
  bits of N, such as the m in question.

  I first parse the rsa encrypted aes key, and given the padding at the beginning
  of the key's plaintext, I can use the ciphertext and partially known plaintext
  to find the unknown bits comprising m.

  I used the coppersmith algorithm as mentioned in Nadia's lecture slide to find m.
  Lecture slides 23 at the top of the bottom page to be exact.

  After retrieving m, I could use AES to decrypt the blocks of the file.
"""

# Gets the encryption of the aes key and the aes encryption of the rest of the file.
def getKeyAndFileEncryption():
  output = open('hw6.pdf.enc.asc', "r").read()
  data = re.search(encrypt_header+"(.*)"+encrypt_footer,output,flags=re.DOTALL).group(1).decode('base64')

  index = 0
  c,  index = parse_mpi(data, index)
  return c, data[index:]


# Determine the known part of the ciphertext, a
# and the ciphertext c
pubkey = RSA.importKey(open('key.pub').read())
k = ceil(pubkey.size()/8)-3-32
print k
message_prefix  = '0001' + 'ff'*k + '00' +'00'*32
a = int(message_prefix, 16)

c, fileEnc = getKeyAndFileEncryption()

# Use coppersmith's method to get a root of the polynomial
# f(m) = (a + m)^3 - c mod N
#      = x^3 + f2*x^2 + f1*x + f0
f2 = 3 * a
f1 = 3 * a**2
f0 = a**3 - c
M = 2**256
N = pubkey.n

# Create matrix as described in Nadia's lecture notes. Slides 23, last page
L = matrix(ZZ, 4, 4, [pow(M, 3), f2*pow(M, 2), f1*M, f0,    # f(x)
                              0,  N*pow(M, 2),    0,  0,    # N*x^2
                              0,            0,  N*M,  0,    # N*x
                              0,            0,    0,  N])   # N
B = L.LLL()
x = PolynomialRing(RationalField(), 'x').gen()

f = sum(B[0][3-i] * x**i / M**i for i in range(4))
roots = f.roots()
m = roots[0][0]
print m

# Having computed m, decrypt the pdf file
aeskey = int_to_binary(m)
iv = Random.new().read(AES.block_size)
cipher = AES.new(aeskey, AES.MODE_CBC, iv)

plaintext = cipher.decrypt(fileEnc) #Is index correct?
f = open("hw6.pdf", "w")
f.write(plaintext)
f.close()
print "Decryption complete!"

