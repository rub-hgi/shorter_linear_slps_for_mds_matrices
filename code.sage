###############################################################################
#### Constructions ############################################################
###############################################################################

def cauchy(a, b):
    a = list(a)
    b = list(b)
    field = a[0].parent()
    rows = []
    for ai in a:
        rows.append([1 / (ai - bi) for bi in b])
    return matrix(field, rows)


def vandermonde(a, b):
    a = list(a)
    b = list(b)
    field = a[0].parent()
    a_rows = []
    for ai in a:
        a_rows.append([ai^i for i in range(len(a))])
    b_rows = []
    for bi in b:
        b_rows.append([bi^i for i in range(len(b))])
    A = matrix(field, a_rows)
    B = matrix(field, b_rows)
    return A * B.inverse()


def hadamard(row):
    row = list(row)
    field = row[0].parent()
    if len(row) == 4:
        a, b, c, d = row
        return matrix(field, 4, 4, [
            [a, b, c, d],
            [b, a, d, c],
            [c, d, a, b],
            [d, c, b, a]])
    elif len(row) == 8:
        a, b, c, d, e, f, g, h = row
        return matrix(field, 8, 8, [
            [a, b, c, d, e, f, g, h],
            [b, a, d, c, f, e, h, g],
            [c, d, a, b, g, h, e, f],
            [d, c, b, a, h, g, f, e],
            [e, f, g, h, a, b, c, d],
            [f, e, h, g, b, a, d, c],
            [g, h, e, f, c, d, a, b],
            [h, g, f, e, d, c, b, a]])


def circ(row):
    row = list(row)
    field = row[0].parent()
    if len(row) == 4:
        a, b, c, d = row
        return matrix(field, 4, 4, [
            [a, b, c, d],
            [d, a, b, c],
            [c, d, a, b],
            [b, c, d, a]])
    if len(row) == 8:
        a, b, c, d, e, f, g, h = row
        return matrix(field, 8, 8, [
            [a, b, c, d, e, f, g, h],
            [h, a, b, c, d, e, f, g],
            [g, h, a, b, c, d, e, f],
            [f, g, h, a, b, c, d, e],
            [e, f, g, h, a, b, c, d],
            [d, e, f, g, h, a, b, c],
            [c, d, e, f, g, h, a, b],
            [b, c, d, e, f, g, h, a]])


def toeplitz(row):
    row = list(row)
    field = row[0].parent()
    if len(row) == 15:
        a, b, c, d, e, f, g, h, t, u, v, w, x, y, z = row
        return matrix(field, 8, 8, [
            [a, b, c, d, e, f, g, h],
            [t, a, b, c, d, e, f, g],
            [u, t, a, b, c, d, e, f],
            [v, u, t, a, b, c, d, e],
            [w, v, u, t, a, b, c, d],
            [x, w, v, u, t, a, b, c],
            [y, x, w, v, u, t, a, b],
            [z, y, x, w, v, u, t, a]])


def subfield(m, sbox_size=0, dim=0):
    """
    Return a binary matrix that multiplies the left and right half of the input with m
    """
    if m[0][0].parent() != GF(2):
        sbox_size = m[0][0].parent().degree()
        dim = m.dimensions()[0]
        binary_m_ij = [[field2bin(m[j][i])
                        for i in range(dim)]
                       for j in range(dim)]
    else:
        assert sbox_size != 0 and dim != 0
        binary_m_ij = [[m.matrix_from_columns(range(j*sbox_size, (j+1)*sbox_size))
                        for j in range(dim)]
                       for m in [m.matrix_from_rows(range(i*sbox_size, (i+1)*sbox_size)) for i in range(dim)]
                      ]

    new_m_ij = [[block_diagonal_matrix(binary_m_ij[j][i], binary_m_ij[j][i])
                 for i in range(dim)]
                for j in range(dim)]

    return block_matrix(GF(2), dim, dim, new_m_ij)


###############################################################################
#### Conversions ##############################################################
###############################################################################

def field2bin(ff_element):
    """
    Converts a finite field element to the binary matrix carrying out the according multiplication
    """
    F = ff_element.parent()
    R = companion_matrix(F.modulus(), format='right')
    L = Integer(ff_element.integer_representation()).digits(base=2, padto=F.degree())
    X = 0
    T = identity_matrix(GF(2), F.degree())
    for i in L:
        if i == 1:
            X = X + T
        T = T * R
    return X


def ff_matrix2bin(ff_matrix):
    """
    Converts a matrix over a finite field to the binary matrix carrying out the according multiplication
    """
    L = []
    for i in ff_matrix:
        for j in i:
            L.append(field2bin(j))
    return block_matrix(ff_matrix.nrows(), ff_matrix.ncols(), L)


def matrix2columns(m, to_file=None):
    """
    converts a matrix to column notation, this is needed by the c implementation
    of the paar algorithms.
    the function can also write this matrix to a correct formatted header file
    that needs to be included in the c implementation
    """
    dim = m[0,0].parent().degree()*m.ncols()
    sbox_size = m[0,0].parent().degree()
    if m[0][0].parent() != GF(2):
        m = ff_matrix2bin(m)

    columns = []
    for col in m.columns():
        columns.append(reduce(lambda a, x: Integer(a)*2 + Integer(x), col))

    if not to_file is None:
        with open(to_file, "w") as f:
            f.write("#define DIM %d\n" % (dim))
            f.write("#define SBOX_SIZE %d\n" % (sbox_size))
            f.write("#define NUM_SBOXES (DIM/SBOX_SIZE)\n\n")
            f.write("uint64_t mds[] = {\n")
            for col in columns[:-1]:
                f.write("    0x%08x,\n" % col)
            f.write("    0x%08x\n" % columns[-1])
            f.write("};")

    return columns


def columns2matrix(cols, dim=16):
    columns = map(lambda x: Integer(list(x)).digits(base=2, padto=dim)[::-1], cols)
    return matrix(GF(2), columns).transpose()


def header2matrix(ifile):
    """
    Reads the matrix stored in the header with name 'ifile' and stores the binary
    matrix in the format needed for the slp heuristics program in the file 'ofile'
    """
    header = ""
    with open(ifile, "r") as f:
        header = map(lambda x: int(x.strip().strip(','), 16), f.readlines()[5:-1])
    binary_matrix = columns2matrix(header, dim=len(header))
    return binary_matrix


def matrix2slp_format(m, to_file=None):
    """
    takes a binary matrix and converts it to the Boyar-Peralta format
    """
    if m[0][0].parent() != GF(2):
        m = ff_matrix2bin(m)

    m = matrix(GF(2), m.dimensions()[0], m.dimensions()[1], m.rows())
    output = "1\n%d %d\n" % m.dimensions()
    for row in m.str().split("\n"):
        output += row[1:-1] + "\n"

    if not to_file is None:
        with open(to_file, "w") as f:
            f.write(output)
    return output


###############################################################################
#### Optimizations ############################################################
###############################################################################

def naive_xor_count(matrix):
    if matrix[0][0].parent() != GF(2):
        matrix = ff_matrix2bin(matrix)
    n, m = matrix.dimensions()
    return matrix.density() * n * m - n


def paar1(m):
    if m[0][0].parent() != GF(2):
        m = ff_matrix2bin(m)

    dim = m.dimensions()[0]
    xor_cnt = naive_xor_count(m)
    m_columns = matrix2columns(m)

    while True:
        hw_max = 0

        for i, col_i in enumerate(m_columns):
            for j, col_j in enumerate(m_columns[i+1:], i+1):
                hw = Integer(col_i & col_j).popcount()
                if hw > hw_max:
                    hw_max = hw
                    i_max = i
                    j_max = j

        if hw_max > 1:
            new_column = m_columns[i_max] & m_columns[j_max]
            m_columns.append(new_column)
            m_columns[i_max] = (new_column^^((1 << dim)-1)) & m_columns[i_max]
            m_columns[j_max] = (new_column^^((1 << dim)-1)) & m_columns[j_max]
            xor_cnt -= hw_max - 1

        if hw_max <= 1:
            break

    return xor_cnt


def slp_heuristic(m, path_to_slp_heuristic_binary="./slp_heuristic", debug=False):
    """
    calls the slp heuristic binary and pipes the matrix in the correct form to
    the cpp program
    """
    from subprocess import Popen, PIPE
    from re import search

    proc = Popen([path_to_slp_heuristic_binary], stdout=PIPE, stdin=PIPE)
    slp_output, _ = proc.communicate(matrix2slp_format(m))

    if debug:
        print(slp_output)

    return int(search(r'\d+',slp_output.split('\n')[0]).group())


def linopt(m, path_to_linopt_binary="./LinOpt", debug=False, nr_times=20000, srand_seed=0):
    """
    calls the LinOpt binary with the correct arguments
    """
    from subprocess import Popen, PIPE, STDOUT
    from tempfile import NamedTemporaryFile
    from re import search

    f = NamedTemporaryFile('w+')
    f.write(matrix2slp_format(m))
    f.seek(0)

    proc = Popen([path_to_linopt_binary, f.name, str(nr_times), str(srand_seed)], stdout=PIPE, stderr=STDOUT)
    slp_output, _ = proc.communicate()
    f.close()

    if debug:
        print(slp_output)

    counts = [line for line in slp_output.split(b'\n') if b'count' in line]
    counts = [int(search(rb'\d+',line).group()) for line in counts]
    counts.sort()

    return counts[0]


###############################################################################
#### Matrices from the Literature #############################################
###############################################################################

FSE_SKOP15_4x4_4_field = GF(2^4, name="a", modulus=PolynomialRing(GF(2), name="x")("x^4 + x + 1"))
FSE_SKOP15_4x4_8_field = GF(2^8, name="a", modulus=PolynomialRing(GF(2), name="x")("x^8 + x^6 + x^5 + x^2 + 1"))
FSE_SKOP15_8x8_4_field = FSE_SKOP15_4x4_4_field
FSE_SKOP15_8x8_8_field = GF(2^8, name="a", modulus=PolynomialRing(GF(2), name="x")("x^8 + x^7 + x^6 + x + 1"))
FSE_SKOP15_4x4_4 = hadamard(map(FSE_SKOP15_4x4_4_field.fetch_int, [0x1, 0x2, 0x8, 0x9]))
FSE_SKOP15_4x4_4_i = hadamard(map(FSE_SKOP15_4x4_4_field.fetch_int, [0x1, 0x4, 0x9, 0xd]))
FSE_SKOP15_4x4_8 = subfield(hadamard(map(FSE_SKOP15_4x4_4_field.fetch_int, [0x1, 0x2, 0x8, 0x9])), 4, 4)
FSE_SKOP15_4x4_8_i = subfield(hadamard(map(FSE_SKOP15_4x4_4_field.fetch_int, [0x1, 0x4, 0x9, 0xd])), 4, 4)
FSE_SKOP15_8x8_4 = hadamard(map(FSE_SKOP15_8x8_4_field.fetch_int, [0x1, 0x2, 0x6, 0x8, 0x9, 0xc, 0xd, 0xa]))
FSE_SKOP15_8x8_4_i = hadamard(map(FSE_SKOP15_8x8_4_field.fetch_int, [0x2, 0x3, 0x4, 0xc, 0x5, 0xa, 0x8, 0xf]))
FSE_SKOP15_8x8_8 = hadamard(map(FSE_SKOP15_8x8_8_field.fetch_int, [0x01, 0x02, 0x03, 0x08, 0x04, 0x91, 0xe1, 0xa9]))
FSE_SKOP15_8x8_8_i = hadamard(map(FSE_SKOP15_8x8_8_field.fetch_int, [0x01, 0x02, 0x03, 0x91, 0x04, 0x70, 0x05, 0xe1]))

C_BeiKraLea16_4x4_4_field = GF(2^4, name="a", modulus=PolynomialRing(GF(2), name="x")("x^4 + x + 1"))
C_BeiKraLea16_4x4_8_field = GF(2^8, name="a", modulus=PolynomialRing(GF(2), name="x")("x^8 + x^6 + x^5 + x + 1"))
C_BeiKraLea16_8x8_8_field = GF(2^8, name="a", modulus=PolynomialRing(GF(2), name="x")("x^8 + x^6 + x^5 + x^2 + 1"))
C_BeiKraLea16_4x4_4 = circ(map(C_BeiKraLea16_4x4_4_field, ["1", "1", "a", "a^-2"]))
C_BeiKraLea16_4x4_8 = circ(map(C_BeiKraLea16_4x4_8_field, ["1", "1", "a", "a^-2"]))
C_BeiKraLea16_8x8_8 = circ(map(C_BeiKraLea16_8x8_8_field, ["1", "1", "a^-1", "a", "a^-1", "a^3", "a^4", "a^-3"]))

FSE_LiuSim16_4x4_4_field = GF(2^4, name="a", modulus=PolynomialRing(GF(2), name="x")("x^4 + x + 1"))
FSE_LiuSim16_4x4_8_field = GF(2^8, name="a", modulus=PolynomialRing(GF(2), name="x")("x^8 + x^7 + x^6 + x + 1"))
FSE_LiuSim16_8x8_8_field = FSE_LiuSim16_4x4_8_field
FSE_LiuSim16_4x4_4 = circ(map(FSE_LiuSim16_4x4_4_field.fetch_int, [0x1, 0x1, 0x9, 0x4]))
FSE_LiuSim16_4x4_8 = circ(map(FSE_LiuSim16_4x4_8_field.fetch_int, [0x01, 0x01, 0x02, 0x91]))
FSE_LiuSim16_8x8_8 = circ(map(FSE_LiuSim16_8x8_8_field.fetch_int, [0x01, 0x01, 0x02, 0xe1, 0x08, 0xe0, 0x01, 0xa9]))

FSE_LiWang16_1 = identity_matrix(GF(2), 4)
FSE_LiWang16_A = matrix(GF(2), [[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,0,0,1]])
FSE_LiWang16_B = matrix(GF(2), [[0,1,1,0],[0,0,1,1],[1,0,0,0],[0,1,0,0]])
FSE_LiWang16_4x4_4 = block_matrix(GF(2), 4, 4, [FSE_LiWang16_1, FSE_LiWang16_1, FSE_LiWang16_A, FSE_LiWang16_B
                                               ,FSE_LiWang16_B, FSE_LiWang16_1, FSE_LiWang16_1, FSE_LiWang16_A
                                               ,FSE_LiWang16_A, FSE_LiWang16_B, FSE_LiWang16_1, FSE_LiWang16_1
                                               ,FSE_LiWang16_1, FSE_LiWang16_A, FSE_LiWang16_B, FSE_LiWang16_1])
FSE_LiWang16_C = matrix(GF(2), [[0,1,1,0],[0,0,0,1],[0,1,0,0],[1,0,0,0]])
FSE_LiWang16_D = FSE_LiWang16_C.inverse() * FSE_LiWang16_C.inverse()
FSE_LiWang16_4x4_4_2 = block_matrix(GF(2), 4, 4, [FSE_LiWang16_C, FSE_LiWang16_1, FSE_LiWang16_1, FSE_LiWang16_1
                                                 ,FSE_LiWang16_1, FSE_LiWang16_1, FSE_LiWang16_D, FSE_LiWang16_C
                                                 ,FSE_LiWang16_1, FSE_LiWang16_C, FSE_LiWang16_1, FSE_LiWang16_D
                                                 ,FSE_LiWang16_1, FSE_LiWang16_D, FSE_LiWang16_C, FSE_LiWang16_1])
FSE_LiWang16_A_i = matrix(GF(2), [[0,1,0,0],[1,0,1,0],[0,0,0,1],[0,1,1,0]])
FSE_LiWang16_B_i = FSE_LiWang16_A_i.inverse()
FSE_LiWang16_C_i = FSE_LiWang16_A_i + FSE_LiWang16_A_i.inverse()
FSE_LiWang16_4x4_4_i = block_matrix(GF(2), 4, 4, [FSE_LiWang16_1  , FSE_LiWang16_A_i, FSE_LiWang16_B_i, FSE_LiWang16_C_i
                                                 ,FSE_LiWang16_A_i, FSE_LiWang16_1  , FSE_LiWang16_C_i, FSE_LiWang16_B_i
                                                 ,FSE_LiWang16_B_i, FSE_LiWang16_C_i, FSE_LiWang16_1  , FSE_LiWang16_A_i
                                                 ,FSE_LiWang16_C_i, FSE_LiWang16_B_i, FSE_LiWang16_A_i, FSE_LiWang16_1])

FSE_LiWang16_8_1 = identity_matrix(GF(2), 8)
FSE_LiWang16_8_A = matrix(GF(2), [[0,1,0,0,0,0,0,0]
                                 ,[0,0,1,0,0,0,0,0]
                                 ,[0,0,0,1,0,0,0,0]
                                 ,[0,0,0,0,1,0,0,0]
                                 ,[0,0,0,0,0,1,0,0]
                                 ,[0,0,0,0,0,0,1,0]
                                 ,[0,0,0,0,0,0,0,1]
                                 ,[1,0,1,0,0,0,0,0]])
FSE_LiWang16_8_B = FSE_LiWang16_8_A.inverse() * FSE_LiWang16_8_A.inverse()
FSE_LiWang16_4x4_8 = block_matrix(GF(2), 4, 4, [FSE_LiWang16_8_1, FSE_LiWang16_8_1, FSE_LiWang16_8_A, FSE_LiWang16_8_B
                                               ,FSE_LiWang16_8_B, FSE_LiWang16_8_1, FSE_LiWang16_8_1, FSE_LiWang16_8_A
                                               ,FSE_LiWang16_8_A, FSE_LiWang16_8_B, FSE_LiWang16_8_1, FSE_LiWang16_8_1
                                               ,FSE_LiWang16_8_1, FSE_LiWang16_8_A, FSE_LiWang16_8_B, FSE_LiWang16_8_1])
FSE_LiWang16_8_C = matrix(GF(2), [[0,0,0,1,0,0,0,0]
                                 ,[0,0,0,0,1,0,0,0]
                                 ,[0,0,0,0,0,1,0,0]
                                 ,[0,0,0,0,0,0,0,1]
                                 ,[0,0,1,0,0,0,0,0]
                                 ,[0,0,0,1,0,0,1,0]
                                 ,[1,0,0,0,0,0,0,0]
                                 ,[0,1,0,0,0,0,0,0]])
FSE_LiWang16_8_D = FSE_LiWang16_8_C.inverse() * FSE_LiWang16_8_C.inverse()
FSE_LiWang16_4x4_8_2 = block_matrix(GF(2), 4, 4, [FSE_LiWang16_8_C, FSE_LiWang16_8_1, FSE_LiWang16_8_1, FSE_LiWang16_8_1
                                                 ,FSE_LiWang16_8_1, FSE_LiWang16_8_1, FSE_LiWang16_8_D, FSE_LiWang16_8_C
                                                 ,FSE_LiWang16_8_1, FSE_LiWang16_8_C, FSE_LiWang16_8_1, FSE_LiWang16_8_D
                                                 ,FSE_LiWang16_8_1, FSE_LiWang16_8_D, FSE_LiWang16_8_C, FSE_LiWang16_8_1])
FSE_LiWang16_8_A_i = matrix(GF(2), [[0,1,0,0,0,0,0,0]
                                   ,[0,0,1,0,0,0,0,0]
                                   ,[0,0,0,1,0,0,0,0]
                                   ,[0,0,0,0,1,0,0,0]
                                   ,[0,0,0,0,0,1,0,0]
                                   ,[0,0,0,0,0,0,1,0]
                                   ,[0,0,0,0,0,0,0,1]
                                   ,[1,0,1,0,0,0,0,0]])
FSE_LiWang16_8_B_i = FSE_LiWang16_8_A_i.inverse()
FSE_LiWang16_8_C_i = FSE_LiWang16_8_A_i + FSE_LiWang16_8_A_i.inverse()
FSE_LiWang16_4x4_8_i = block_matrix(GF(2), 4, 4, [FSE_LiWang16_8_1  , FSE_LiWang16_8_A_i, FSE_LiWang16_8_B_i, FSE_LiWang16_8_C_i
                                                 ,FSE_LiWang16_8_A_i, FSE_LiWang16_8_1  , FSE_LiWang16_8_C_i, FSE_LiWang16_8_B_i
                                                 ,FSE_LiWang16_8_B_i, FSE_LiWang16_8_C_i, FSE_LiWang16_8_1  , FSE_LiWang16_8_A_i
                                                 ,FSE_LiWang16_8_C_i, FSE_LiWang16_8_B_i, FSE_LiWang16_8_A_i, FSE_LiWang16_8_1])
FSE_LiWang16_8_A_i_2 = matrix(GF(2), [[1,0,0,0,0,0,0,0]
                                     ,[0,1,0,0,0,0,0,0]
                                     ,[1,0,1,0,0,0,0,0]
                                     ,[1,1,0,1,0,0,0,0]
                                     ,[0,0,0,0,0,1,0,0]
                                     ,[0,0,0,0,1,0,0,0]
                                     ,[0,0,0,0,0,0,0,1]
                                     ,[0,0,0,0,0,0,1,0]])
FSE_LiWang16_8_B_i_2 = matrix(GF(2), [[0,0,0,0,0,0,1,1]
                                     ,[1,0,0,0,0,0,0,0]
                                     ,[0,0,0,0,0,0,1,0]
                                     ,[0,0,1,0,0,0,0,1]
                                     ,[0,1,0,1,0,0,0,0]
                                     ,[1,0,0,1,0,0,0,0]
                                     ,[0,0,0,0,0,1,0,0]
                                     ,[0,0,0,0,1,0,0,0]])
FSE_LiWang16_8_C_i_2 = matrix(GF(2), [[0,0,0,0,1,0,0,0]
                                     ,[0,0,0,0,0,0,0,1]
                                     ,[0,1,0,0,0,1,0,0]
                                     ,[0,0,0,0,0,0,1,0]
                                     ,[1,0,0,0,0,0,0,0]
                                     ,[0,0,1,0,0,0,0,1]
                                     ,[0,0,0,1,0,0,0,0]
                                     ,[0,1,0,0,0,0,0,0]])
FSE_LiWang16_4x4_8_i_2 = block_matrix(GF(2), 4, 4, [FSE_LiWang16_8_1    , FSE_LiWang16_8_A_i_2, FSE_LiWang16_8_B_i_2, FSE_LiWang16_8_C_i_2
                                                   ,FSE_LiWang16_8_C_i_2, FSE_LiWang16_8_1    , FSE_LiWang16_8_A_i_2, FSE_LiWang16_8_B_i_2
                                                   ,FSE_LiWang16_8_B_i_2, FSE_LiWang16_8_C_i_2, FSE_LiWang16_8_1    , FSE_LiWang16_8_A_i_2
                                                   ,FSE_LiWang16_8_A_i_2, FSE_LiWang16_8_B_i_2, FSE_LiWang16_8_C_i_2, FSE_LiWang16_8_1])

ToSC_SarSye16_4x4_4_field = GF(2^4, name="a", modulus=PolynomialRing(GF(2), name="x")("x^4 + x^3 + 1"))
ToSC_SarSye16_4x4_8_field = GF(2^8, name="a", modulus=PolynomialRing(GF(2), name="x")("x^8 + x^7 + x^6 + x + 1"))
ToSC_SarSye16_4x4_4 = matrix(ToSC_SarSye16_4x4_4_field, map(lambda x: map(ToSC_SarSye16_4x4_4_field, x), [["1", "1", "a", "a^-1"], ["a^-2", "1", "1", "a"], ["1", "a^-2", "1", "1"], ["a^-1", "1", "a^-2", "1"]]))
ToSC_SarSye16_4x4_8 = matrix(ToSC_SarSye16_4x4_8_field, map(lambda x: map(ToSC_SarSye16_4x4_8_field, x), [["1", "1", "a", "a^-1"], ["a^-2", "1", "1", "a"], ["1", "a^-2", "1", "1"], ["a^-1", "1", "a^-2", "1"]]))

ToSC_SarSye16_4x4_4_i_field = GF(2^4, name="a", modulus=PolynomialRing(GF(2), name="x")("x^4 + x + 1"))
ToSC_SarSye16_4x4_8_i_field = GF(2^8, name="a", modulus=PolynomialRing(GF(2), name="x")("x^8 + x^6 + x^5 + x^2 + 1"))
ToSC_SarSye16_4x4_4_i = matrix(ToSC_SarSye16_4x4_4_i_field, map(lambda x: map(ToSC_SarSye16_4x4_4_i_field, x), [["1", "a", "a^2", "1"], ["a", "1", "1", "a^2"], ["a^3", "a", "1", "a"], ["a", "a^3", "a", "1"]]))
ToSC_SarSye16_4x4_8_i = matrix(ToSC_SarSye16_4x4_8_i_field, map(lambda x: map(ToSC_SarSye16_4x4_8_i_field, x), [["1", "a", "1", "a^211"], ["a", "1", "a^211", "1"], ["a^-2", "a^209", "1", "a"], ["a^209", "a^-2", "a", "1"]]))

ACISP_SarSye17_8x8_4_field = GF(2^4, name="a", modulus=PolynomialRing(GF(2), name="x")("x^4 + x + 1"))
ACISP_SarSye17_8x8_8_field = GF(2^8, name="a", modulus=PolynomialRing(GF(2), name="x")("x^8 + x^7 + x^6 + x + 1"))
ACISP_SarSye17_8x8_4 = toeplitz(map(ACISP_SarSye17_8x8_4_field, ["a", "1", "a^4", "1", "a^5", "a^14", "a^7", "a^8", "a^3", "a^6", "a^14", "a^14", "a^8", "a^6", "a^3"]))
ACISP_SarSye17_8x8_8 = toeplitz(map(ACISP_SarSye17_8x8_8_field, ["1", "1", "a", "a^253", "1", "a^253", "a^252", "a^157", "a^158", "a^253", "a^254", "a", "a^254", "a^2", "a"]))

ePrint_JeaPeySim_field = GF(2^4, name="a", modulus=PolynomialRing(GF(2), name="x")("x^4 + x + 1"))
ePrint_JeaPeySim_field_8 = GF(2^8, name="a", modulus=PolynomialRing(GF(2), name="x")("x^8 + x^7 + x^6 + x + 1"))
ePrint_JeaPeySim_4x4_4   = matrix(ePrint_JeaPeySim_field, map(lambda x: map(ePrint_JeaPeySim_field.fetch_int, x), [[0x1, 0x1, 0x1, 0x2], [0x1, 0x2, 0xd, 0x1], [0x2, 0xd, 0x1, 0x1], [0xd, 0x1, 0x2, 0x1]]))
ePrint_JeaPeySim_4x4_4_i = matrix(ePrint_JeaPeySim_field, map(lambda x: map(ePrint_JeaPeySim_field.fetch_int, x), [[0x2, 0x1, 0x1, 0x9], [0x1, 0x4, 0xf, 0x1], [0xd, 0x9, 0x4, 0x1], [0x1, 0xd, 0x1, 0x2]]))
ePrint_JeaPeySim_4x4_8   = subfield(ePrint_JeaPeySim_4x4_4)
ePrint_JeaPeySim_4x4_8_i = subfield(ePrint_JeaPeySim_4x4_4_i)
ePrint_JeaPeySim_8x8_4   = hadamard(map(ePrint_JeaPeySim_field.fetch_int, [0x1, 0x2, 0x6, 0x8, 0x9, 0xc, 0xd, 0xa]))
ePrint_JeaPeySim_8x8_4_i = hadamard(map(ePrint_JeaPeySim_field.fetch_int, [0x2, 0x3, 0x4, 0xc, 0x5, 0xa, 0x8, 0xf]))
ePrint_JeaPeySim_8x8_8   = circ(map(ePrint_JeaPeySim_field_8.fetch_int, [0x1, 0x1, 0x2, 0xe1, 0x8, 0xe0, 0x1, 0xa9]))
ePrint_JeaPeySim_8x8_8_i = hadamard(map(ePrint_JeaPeySim_field_8.fetch_int, [0x1, 0x2, 0x4, 0x91, 0x6a, 0x5b, 0xe1, 0xa9]))

F_13 = GF(2^4, name="a", modulus=PolynomialRing(GF(2), name="x")("x^4 + x + 1"))
F_11d = GF(2^8, name="a", modulus=PolynomialRing(GF(2), name="x")("x^8 + x^4 + x^3 + x^2 + 1"))
F_11b = GF(2^8, name="a", modulus=PolynomialRing(GF(2), name="x")("x^8 + x^4 + x^3 + x + 1"))
F_169 = GF(2^8, name="a", modulus=PolynomialRing(GF(2), name="x")("x^8 + x^6 + x^5 + x^3 + 1"))

Twofish = matrix(F_169, 4, 4, map(lambda x: map(F_169.fetch_int, x), [
    [0x01, 0xef, 0x5b, 0x5b],
    [0x5b, 0xef, 0xef, 0x01],
    [0xef, 0x5b, 0x01, 0xef],
    [0xef, 0x01, 0xef, 0x5b]]))
Clefia_M0 = hadamard(map(F_11d.fetch_int, [1, 2, 4, 6]))
Clefia_M1 = hadamard(map(F_11d.fetch_int, [1, 8, 2, 10]))

Whirlwind_M0 = hadamard(map(F_13.fetch_int, [0x5, 0x4, 0xa, 0x6, 0x2, 0xd, 0x8, 0x3]))
Whirlwind_M1 = hadamard(map(F_13.fetch_int, [0x5, 0xe, 0x4, 0x7, 0x1, 0x3, 0xf, 0x8]))
Whirlpool = circ(map(F_11d.fetch_int, [0x1, 0x1, 0x4, 0x1, 0x8, 0x5, 0x2, 0x9]))
Khazad = hadamard(map(F_11d.fetch_int, [0x1, 0x3, 0x4, 0x5, 0x6, 0x8, 0xb, 0x7]))
Grostl = circ(map(F_11b.fetch_int, [0x2, 0x2, 0x3, 0x4, 0x5, 0x3, 0x5, 0x7]))
Joltik = hadamard(map(F_13.fetch_int, [1, 4, 9, 13]))
Anubis = hadamard(map(F_11d.fetch_int, [1, 2, 4, 6]))
Photon = matrix(F_11b, 4, 4, [[0,1,0,0],[0,0,1,0],[0,0,0,1],map(F_11b.fetch_int, [1,2,1,4])])
Fox_Mu4 = matrix(F_11b, 4, 4, map(lambda x: map(F_11b.fetch_int, x), [[1,1,1,2], [1, 0xfe, 2, 1], [0xfe, 2, 1, 1], [2, 1, 0xfe, 1]]))
Fox_Mu8 = matrix(F_11b, 8, 8, map(lambda x: map(F_11b.fetch_int, x), [
    [0x01,0x01,0x01,0x01,0x01,0x01,0x01,0x03],
    [0x01,0x03,0x82,0x02,0x04,0xfc,0x7e,0x01],
    [0x03,0x82,0x02,0x04,0xfc,0x7e,0x01,0x01],
    [0x82,0x02,0x04,0xfc,0x7e,0x01,0x03,0x01],
    [0x02,0x04,0xfc,0x7e,0x01,0x03,0x82,0x01],
    [0x04,0xfc,0x7e,0x01,0x03,0x82,0x02,0x01],
    [0xfc,0x7e,0x01,0x03,0x82,0x02,0x04,0x01],
    [0x7e,0x01,0x03,0x82,0x02,0x04,0xfc,0x01]]))

AES_field = GF(2^8, name="a", modulus=PolynomialRing(GF(2), name="x")("x^8 + x^4 + x^3 + x + 1"))
AES = circ(map(AES_field.fetch_int, [0x2, 0x3, 0x1, 0x1]))
SmallScale_AES = circ(map(F_13.fetch_int, [0x2, 0x3, 0x1, 0x1]))

def ARIA():
    def ll(x):
        x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15 = x
        y0  = x3 + x4 + x6 + x8  + x9  + x13 + x14
        y1  = x2 + x5 + x7 + x8  + x9  + x12 + x15
        y2  = x1 + x4 + x6 + x10 + x11 + x12 + x15
        y3  = x0 + x5 + x7 + x10 + x11 + x13 + x14
        y4  = x0 + x2 + x5 + x8  + x11 + x14 + x15
        y5  = x1 + x3 + x4 + x9  + x10 + x14 + x15
        y6  = x0 + x2 + x7 + x9  + x10 + x12 + x13
        y7  = x1 + x3 + x6 + x8  + x11 + x12 + x13
        y8  = x0 + x1 + x4 + x7  + x10 + x13 + x15
        y9  = x0 + x1 + x5 + x6  + x11 + x12 + x14
        y10 = x2 + x3 + x5 + x6  + x8  + x13 + x15
        y11 = x2 + x3 + x4 + x7  + x9  + x12 + x14
        y12 = x1 + x2 + x6 + x7  + x9  + x11 + x12
        y13 = x0 + x3 + x6 + x7  + x8  + x10 + x13
        y14 = x0 + x3 + x4 + x5  + x9  + x11 + x14
        y15 = x1 + x2 + x4 + x5  + x8  + x10 + x15
        return VectorSpace(GF(2), 16)([y0,y1,y2,y3,y4,y5,y6,y7,y8,y9,y10,y11,y12,y13,y14,y15])
    return matrix(GF(2^8), [ll(Integer(1 << i).digits(base=2, padto=16)) for i in range(16)])
ARIA = ARIA()

QARMA64_nul = matrix(GF(2), [[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])
QARMA64_rho = matrix(GF(2), [[0,0,0,1],[1,0,0,0],[0,1,0,0],[0,0,1,0]])
QARMA64_rhosq = QARMA64_rho * QARMA64_rho
QARMA64 = block_matrix(GF(2), 4, 4,
    [QARMA64_nul, QARMA64_rho, QARMA64_rhosq, QARMA64_rho,
     QARMA64_rho, QARMA64_nul, QARMA64_rho, QARMA64_rhosq,
     QARMA64_rhosq, QARMA64_rho, QARMA64_nul, QARMA64_rho,
     QARMA64_rho, QARMA64_rhosq, QARMA64_rho, QARMA64_nul])

QARMA128_nul = matrix(GF(2), [[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0],[0,0,0,0,0,0,0,0]])
QARMA128_rho = matrix(GF(2), [[0,0,0,0,0,0,0,1],[1,0,0,0,0,0,0,0],[0,1,0,0,0,0,0,0],[0,0,1,0,0,0,0,0],[0,0,0,1,0,0,0,0],[0,0,0,0,1,0,0,0],[0,0,0,0,0,1,0,0],[0,0,0,0,0,0,1,0]])
QARMA128_rho4 = QARMA128_rho * QARMA128_rho * QARMA128_rho * QARMA128_rho
QARMA128_rho5 = QARMA128_rho4 * QARMA128_rho
QARMA128 = block_matrix(GF(2), 4, 4,
    [QARMA128_nul, QARMA128_rho, QARMA128_rho4, QARMA128_rho5,
     QARMA128_rho5, QARMA128_nul, QARMA128_rho, QARMA128_rho4,
     QARMA128_rho4, QARMA128_rho5, QARMA128_nul, QARMA128_rho,
     QARMA128_rho, QARMA128_rho4, QARMA128_rho5, QARMA128_nul])

MIDORI = circ(map(F_13.fetch_int, [0, 1, 1, 1]))

SKINNY = matrix(F_13, map(lambda x: map(F_13.fetch_int, x), [[0x1, 0x0, 0x1, 0x1], [0x1, 0x0, 0x0, 0x0], [0x0, 0x1, 0x1, 0x0], [0x1, 0x0, 0x1, 0x0]]))

PRINCE_M0 = matrix(GF(2), [[0,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
PRINCE_M1 = matrix(GF(2), [[1,0,0,0],[0,0,0,0],[0,0,1,0],[0,0,0,1]])
PRINCE_M2 = matrix(GF(2), [[1,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,0,1]])
PRINCE_M3 = matrix(GF(2), [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,0]])
PRINCE_M_0 = block_matrix(GF(2), 4, 4,
                     [PRINCE_M0, PRINCE_M1, PRINCE_M2, PRINCE_M3,
                      PRINCE_M1, PRINCE_M2, PRINCE_M3, PRINCE_M0,
                      PRINCE_M2, PRINCE_M3, PRINCE_M0, PRINCE_M1,
                      PRINCE_M3, PRINCE_M0, PRINCE_M1, PRINCE_M2])
PRINCE_M_1 = block_matrix(GF(2), 4, 4,
                     [PRINCE_M1, PRINCE_M2, PRINCE_M3, PRINCE_M0,
                      PRINCE_M2, PRINCE_M3, PRINCE_M0, PRINCE_M1,
                      PRINCE_M3, PRINCE_M0, PRINCE_M1, PRINCE_M2,
                      PRINCE_M0, PRINCE_M1, PRINCE_M2, PRINCE_M3])

PRIDE_L_0 = matrix(GF(2),
    [[0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0],
     [0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0],
     [0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0],
     [0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1],
     [1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0],
     [0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0],
     [0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0],
     [0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1],
     [1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0],
     [0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0],
     [0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0],
     [0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1],
     [1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0],
     [0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0],
     [0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0],
     [0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0]])

PRIDE_L_1 = matrix(GF(2),
    [[1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
     [0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0],
     [0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0],
     [0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0],
     [0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1],
     [0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0],
     [1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0],
     [1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0],
     [0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0],
     [0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0],
     [0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1],
     [0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1],
     [0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0],
     [0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0],
     [0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0]])

PRIDE_L_2 = matrix(GF(2),
    [[0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,1],
     [0,0,0,0,0,1,1,0,1,0,0,0,0,0,0,0],
     [0,0,0,0,0,0,1,1,0,1,0,0,0,0,0,0],
     [1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0],
     [1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0],
     [0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0],
     [0,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0],
     [0,0,0,1,1,0,0,0,0,0,0,0,0,0,1,0],
     [0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1],
     [0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0],
     [0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0],
     [0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0],
     [1,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0],
     [0,1,0,0,0,0,0,0,0,0,0,0,1,1,0,0],
     [0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,0],
     [0,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1]])

PRIDE_L_3 = matrix(GF(2),
    [[1,0,0,0,1,0,0,0,0,0,0,0,1,0,0,0],
     [0,1,0,0,0,1,0,0,0,0,0,0,0,1,0,0],
     [0,0,1,0,0,0,1,0,0,0,0,0,0,0,1,0],
     [0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,1],
     [1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0],
     [0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0],
     [0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0],
     [0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0],
     [0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0],
     [0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0],
     [0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0],
     [0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1],
     [1,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0],
     [0,1,0,0,0,0,0,0,0,1,0,0,0,1,0,0],
     [0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0],
     [0,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1]])


###############################################################################
#### Our Matrices #############################################################
###############################################################################

M_4_4_1 = matrix(GF(2), [[0,0,0,1], [0,0,1,0], [0,1,0,0], [1,0,0,0]])
M_4_4_2 = matrix(GF(2), [[0,0,1,1], [1,0,0,1], [1,1,0,0], [0,1,0,0]])
M_4_4_3 = matrix(GF(2), [[1,1,0,1], [1,1,0,0], [0,1,0,1], [0,0,1,0]])
M_4_4_4 = matrix(GF(2), [[1,1,0,0], [0,1,0,1], [1,0,1,1], [0,0,0,1]])
M_4_4 = block_matrix(GF(2), 4, 4, [M_4_4_1, M_4_4_2, M_4_4_3, M_4_4_4
                            ,M_4_4_2, M_4_4_1, M_4_4_4, M_4_4_3
                            ,M_4_4_3, M_4_4_4, M_4_4_1, M_4_4_2
                            ,M_4_4_4, M_4_4_3, M_4_4_2, M_4_4_1])

M_4_8 = subfield(M_4_4, 4, 4)
M_8_4 = cauchy(map(F_13, ["a^3 + a^2", "a^3", "a^3 + a + 1", "a + 1", "0", "a^3 + a^2 + a + 1", "a^2", "a^2 + a + 1"])
              ,map(F_13, ["1", "a^2 + 1", "a^3 + a^2 + a", "a^3 + 1", "a^3 + a^2 + 1", "a^2 + a", "a^3 + a", "a"]))
M_8_8 = subfield(M_8_4, 4, 8)
M_i_4_8 = subfield(ToSC_SarSye16_4x4_4_i, 4, 4)
M_i_8_4 = vandermonde(map(F_13, ["a^3 + a + 1", "a + 1", "a^3 + a^2 + a", "a^3 + a^2 + 1", "a^3 + 1", "a^3", "0", "a^3 + a"])
                     ,map(F_13, ["a^2 + a + 1", "a^3 + a^2 + a + 1", "a", "1", "a^2 + 1", "a^2", "a^3 + a^2", "a^2 + a"]))
M_i_8_8 = subfield(M_i_8_4, 4, 8)


###############################################################################
#### Matrix Collections #######################################################
###############################################################################

matrices_4x4_4 = {
    "FSE_SKOP15_4x4_4" : FSE_SKOP15_4x4_4,
    "FSE_LiuSim16_4x4_4" : FSE_LiuSim16_4x4_4,
    "FSE_LiWang16_4x4_4" : FSE_LiWang16_4x4_4,
    "FSE_LiWang16_4x4_4_2" : FSE_LiWang16_4x4_4_2,
    "C_BeiKraLea16_4x4_4" : C_BeiKraLea16_4x4_4,
    "ToSC_SarSye16_4x4_4" : ToSC_SarSye16_4x4_4,
    "ePrint_JeaPeySim_4x4_4" : ePrint_JeaPeySim_4x4_4,
    "SmallScale_AES" : SmallScale_AES,
    "M_4_4" : M_4_4}

matrices_i_4x4_4 = {
    "FSE_SKOP15_i_4x4_4" : FSE_SKOP15_4x4_4_i,
    "FSE_LiWang16_i_4x4_4" : FSE_LiWang16_4x4_4_i,
    "ToSC_SarSye16_i_4x4_4" : ToSC_SarSye16_4x4_4_i,
    "ePrint_JeaPeySim_i_4x4_4" : ePrint_JeaPeySim_4x4_4_i,
    "Joltik" : Joltik}

matrices_4x4_8 = {
    "FSE_SKOP15_4x4_8" : FSE_SKOP15_4x4_8,
    "FSE_LiuSim16_4x4_8" : FSE_LiuSim16_4x4_8,
    "FSE_LiWang16_4x4_8" : FSE_LiWang16_4x4_8,
    "FSE_LiWang16_4x4_8_2" : FSE_LiWang16_4x4_8_2,
    "C_BeiKraLea16_4x4_8" : C_BeiKraLea16_4x4_8,
    "ToSC_SarSye16_4x4_8" : ToSC_SarSye16_4x4_8,
    "ePrint_JeaPeySim_4x4_8" : ePrint_JeaPeySim_4x4_8,
    "AES" : AES,
    "Clefia_M0" : Clefia_M0,
    "Clefia_M1" : Clefia_M1,
    "Fox_Mu4" : Fox_Mu4,
    "Twofish" : Twofish,
    "M_4_8" : M_4_8}

matrices_i_4x4_8 = {
    "FSE_SKOP15_i_4x4_8" : FSE_SKOP15_4x4_8_i,
    "FSE_LiWang16_i_4x4_8" : FSE_LiWang16_4x4_8_i,
    "FSE_LiWang16_i_4x4_8_2" : FSE_LiWang16_4x4_8_i_2,
    "ToSC_SarSye16_i_4x4_8" : ToSC_SarSye16_4x4_8_i,
    "ePrint_JeaPeySim_i_4x4_8" : ePrint_JeaPeySim_4x4_8_i,
    "Anubis" : Anubis,
    "M_i_4_8" : M_i_4_8}

matrices_8x8_4 = {
    "FSE_SKOP15_8x8_4" : FSE_SKOP15_8x8_4,
    "ACISP_SarSye17_8x8_4" : ACISP_SarSye17_8x8_4,
    "Whirlwind_M0" : Whirlwind_M0,
    "Whirlwind_M1" : Whirlwind_M1,
    "M_8_4" : M_8_4}

matrices_i_8x8_4 = {
    "FSE_SKOP15_i_8x8_4" : FSE_SKOP15_8x8_4_i,
    "M_i_8_4" : M_i_8_4}

matrices_8x8_8 = {
    "FSE_SKOP15_8x8_8" : FSE_SKOP15_8x8_8,
    "FSE_LiuSim16_8x8_8" : FSE_LiuSim16_8x8_8,
    "C_BeiKraLea16_8x8_8" : C_BeiKraLea16_8x8_8,
    "ACISP_SarSye17_8x8_8" : ACISP_SarSye17_8x8_8,
    "Fox_Mu8" : Fox_Mu8,
    "Grostl" : Grostl,
    "Whirlpool" : Whirlpool,
    "M_8x8_8" : M_8_8,
}

matrices_i_8x8_8 = {
    "FSE_SKOP15_i_8x8_8" : FSE_SKOP15_8x8_8_i,
    "ePrint_JeaPeySim_i_8x8_8" : ePrint_JeaPeySim_4x4_8_i,
    "Khazad" : Khazad,
    "M_i_8x8_8" : M_i_8_8,
}

matrices_non_MDS = {
    "QARMA64" : QARMA64,
    "QARMA128" : QARMA128,
    "MIDORI" : MIDORI,
    "SKINNY" : SKINNY,
    "PRINCE_M_0" : PRINCE_M_0,
    "PRINCE_M_1" : PRINCE_M_1,
    "PRIDE_L_0" : PRIDE_L_0,
    "PRIDE_L_1" : PRIDE_L_1,
    "PRIDE_L_2" : PRIDE_L_2,
    "PRIDE_L_3" : PRIDE_L_3,
}

all_matrices = matrices_4x4_4.copy()
all_matrices.update(matrices_i_4x4_4)
all_matrices.update(matrices_4x4_8)
all_matrices.update(matrices_i_4x4_8)
all_matrices.update(matrices_8x8_4)
all_matrices.update(matrices_i_8x8_4)
all_matrices.update(matrices_8x8_8)
all_matrices.update(matrices_i_8x8_8)
all_matrices.update(matrices_non_MDS)
