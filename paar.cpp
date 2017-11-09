#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <vector>

using namespace std;

void prepare_matrix(uint64_t *augmented_matrix, uint64_t *matrix) {
    for (uint64_t i=0; i<DIM; i++) {
        matrix[i] = matrix[i] & ((1llu<<DIM)-1);
        augmented_matrix[i] = matrix[i];
    }
}

int is_invertible(uint64_t *matrix) {
    // check all 2^dim linear combination
    for (uint64_t coefficients=1; coefficients<=(1llu<<DIM)-1; coefficients++) {
        uint64_t tmp_coefficients;
        uint64_t sum = 0x0;
        tmp_coefficients = coefficients;
        for (uint64_t i=0; i<DIM; i++) {
            if ( (tmp_coefficients & (1llu<<(DIM-1))) != 0 ) {
                sum = sum^matrix[i];
            }
            tmp_coefficients = tmp_coefficients << 1;
        }
        if (sum == 0) {
            return 0;
        }
    }
    return 1;
}

// naive MDS check by computing branch number
int is_MDS(uint64_t *matrix) {
    for (uint64_t input=1; input<=(1llu<<DIM)-1; input++) {
        uint64_t output = 0x0;
        uint64_t tmp_input = input;
        for (int i=0; i<DIM; i++) {
            if ( (tmp_input & (1ll<<(DIM-1))) != 0 ) {
                output = output^matrix[i];
            }
            tmp_input = tmp_input << 1;
        }

        // check if MDS property still holds
        uint64_t branch_number = 0;
        tmp_input = input;
        for (int i=0; i<NUM_SBOXES; i++) {
            if ( (tmp_input & ((1<<SBOX_SIZE)-1)) != 0 ) {
                branch_number++;
            }
            tmp_input = tmp_input >> SBOX_SIZE;
        }
        for (int i=0; i<NUM_SBOXES; i++) {
            if ( (output & ((1<<SBOX_SIZE)-1)) != 0 ) {
                branch_number++;
            }
            output = output >> SBOX_SIZE;
        }
        if (branch_number < NUM_SBOXES+1) {
            return 0;
        }
    }
    return 1;
}

int hamming_weight(uint64_t input) {
    return __builtin_popcount(input & 0xffffffff) + __builtin_popcount((input >> 32) & 0xffffffff);
}

int paar_algorithm1(uint64_t *input_matrix) {
    int xor_count = 0;
    uint64_t number_of_columns = DIM;
    int hw_max;
    int i_max = 0, j_max = 0;
    uint64_t tmp;
    int hw;
    uint64_t new_column;
    vector<pair<int, int>> program;

    // compute naive xor count
    for (uint64_t i=0; i<DIM; i++) {
        xor_count += hamming_weight(input_matrix[i]);
    }
    xor_count -= DIM;
    cout << "Naive XOR count: " << xor_count << endl;
    cout << "SLP:" << endl << endl;

    do {
        hw_max = 0;
        for (uint64_t i=0; i<number_of_columns; i++) {
            for (uint64_t j=i+1; j<number_of_columns; j++) {
                tmp = input_matrix[i] & input_matrix[j];
                hw = hamming_weight(tmp);
                if (hw > hw_max) {
                    hw_max = hw;
                    i_max = i;
                    j_max = j;
                }
            }
        }
        if (hw_max > 1) {
            new_column = input_matrix[i_max] & input_matrix[j_max];
            input_matrix[number_of_columns] = new_column;
            input_matrix[i_max] = ( new_column^((1llu << DIM)-1) ) & input_matrix[i_max];
            input_matrix[j_max] = ( new_column^((1llu << DIM)-1) ) & input_matrix[j_max];
            xor_count -= (hw_max-1);
            number_of_columns++;
            program.push_back(make_pair(i_max, j_max));
        }
    } while (hw_max > 1);

    int ctr = DIM;
    for (auto const& prog:program) {
        cout << "x" << ctr << " = x" << prog.first << " + x" << prog.second << endl;
        ctr++;
    }

    for (uint64_t i=0; i<DIM; i++) {
        bool plus_flag = 0;
        cout << endl << "y" << i;
        for (uint64_t j=0; j<number_of_columns; j++) {
            if ( ( input_matrix[j] & (1ll<<(DIM-1-i)) ) != 0) {
                if (plus_flag == 0) {
                    cout << " = x" << j;
                    plus_flag = 1;
                }
                else {
                    cout << " + x" << j;
                }
            }
        }
    }

    cout << endl << endl;

    return xor_count;
}

int global_best_xor_count;

int paar_algorithm_2_recursive(uint64_t *input_matrix, int current_xor_count, int number_of_columns, vector<pair<int, int>> program) {
    int xor_count_new;
    int xor_count_best = current_xor_count;
    int xor_count_tmp;
    uint64_t tmp;
    int hw_max;
    int hw;
    int i, j;
    uint64_t new_column;

    int i_max_list[1000];
    int j_max_list[1000];
    int list_end = 0;

    uint64_t *input_matrix_new;
    input_matrix_new = (uint64_t *)malloc((DIM+200)*sizeof(uint64_t));

    hw_max = 0;
    for (i=0; i<number_of_columns; i++) {
        for (j=i+1; j<number_of_columns; j++) {
            tmp = input_matrix[i] & input_matrix[j];
            hw = hamming_weight(tmp);
            if (hw > hw_max) {
                hw_max = hw;
                i_max_list[0] = i;
                j_max_list[0] = j;
                list_end = 1;
            }
            else if ((hw == hw_max) && (hw != 1)) {
                i_max_list[list_end] = i;
                j_max_list[list_end] = j;
                list_end++;
            }
        }
    }
    if (hw_max > 1) {
        for (i=0; i<list_end; i++) {
            for (j=0; j<number_of_columns; j++) {
                input_matrix_new[j] = input_matrix[j];
            }

            vector<pair<int, int>> program_new;
            for (auto const& prog:program) {
                program_new.push_back(prog);
            }

            new_column = input_matrix_new[i_max_list[i]] & input_matrix_new[j_max_list[i]];
            input_matrix_new[number_of_columns] = new_column;
            input_matrix_new[i_max_list[i]] = ( new_column^((1llu << DIM)-1) ) & input_matrix_new[i_max_list[i]];
            input_matrix_new[j_max_list[i]] = ( new_column^((1llu << DIM)-1) ) & input_matrix_new[j_max_list[i]];
            xor_count_new = current_xor_count - (hw_max-1);
            program_new.push_back(make_pair(i_max_list[i], j_max_list[i]));
            xor_count_tmp = paar_algorithm_2_recursive(input_matrix_new, xor_count_new, number_of_columns+1, program_new);
            if (xor_count_tmp < xor_count_best) {
                xor_count_best = xor_count_tmp;
            }
        }
    }

    free(input_matrix_new);

    if (xor_count_best < global_best_xor_count) {
        global_best_xor_count = xor_count_best;

        int ctr = DIM;
        for (auto const& prog:program) {
            cout << "x" << ctr << " = x" << prog.first << " + x" << prog.second << endl;
            ctr++;
        }

        for (int i=0; i<DIM; i++) {
            int plus_flag = 0;
            cout << endl << "y" << i;
            for (int j=0; j<number_of_columns; j++) {
                if ( ( input_matrix[j] & (1ll<<(DIM-1-i)) ) != 0) {
                    if (plus_flag == 0) {
                        cout << " = x" << j;
                        plus_flag = 1;
                    }
                    else {
                        cout << " + x" << j;
                    }
                }
            }
        }

        cout << endl << endl << "xor count (tmp) = " << global_best_xor_count << endl << endl << endl;
    }

    return xor_count_best;
}

int paar_algorithm2(uint64_t *input_matrix) {
    int xor_count = 0;

    // compute naive xor count
    for (int i=0; i<DIM; i++) {
        xor_count += hamming_weight(input_matrix[i]);
    }
    xor_count -= DIM;

    global_best_xor_count = xor_count;
    cout << "xor count (start) = " << global_best_xor_count << endl << endl;

    vector<pair<int, int>> program;

    xor_count = paar_algorithm_2_recursive(input_matrix, xor_count, DIM, program);

    return xor_count;
}

int main(int argc, char ** argv){
    // store matrix column-wise
    uint64_t *input_matrix;
    input_matrix = (uint64_t *)malloc((DIM+200)*sizeof(uint64_t));
    prepare_matrix(input_matrix, mds);

#ifdef PAAR1
    cout << "Paar1 XOR count:" << paar_algorithm1(input_matrix) << endl;
#elif PAAR2
    cout << "Paar2 XOR count:" << paar_algorithm2(input_matrix) << endl;
#endif

    free(input_matrix);

    return 0;
}
