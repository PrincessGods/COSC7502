/*
 * ga_island.c by Jordan T. Bishop 2019
 * Developed for COSC3500/7502 at the University of Queensland
 *
 * Compile with:
 * mpicc -std=gnu99 ga_island.c -lm -o ga_island.exe -O3 -Wall -Wextra
 *
 * Run with ga_island.sh Slurm submission script.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <stdarg.h>
#include <mpi.h>

// struct for use with qsort comparator
typedef struct {
    int pop_idx;
    int fitness;
} pop_fitness_map;

const int root = 0;

/* The C "%" operator is actually a remainder operator, so we need our own
 * function to do proper modular arithmetic.
 */
int modulo(int a, int b){
    return (a % b + b) % b;
}

// Print formatted message on given rank
void rank_print(int my_rank, char* format, ...) {
    va_list args;
    va_start(args, format); 
    printf("Rank %d: ", my_rank);
    vprintf(format, args);
    va_end(args);
}

// Print formatted message on root rank only
void root_print(int my_rank, char* format, ...) {
    va_list args;
    if (my_rank != root) {
        return;
    }
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
}

/* Fill the pop array with random zeros/ones
 * Algorithm 21 from 'Essentials of Metaheuristics'
 */
void init_population(int* pop, int pop_size, int prob_size) {
    for (int i = 0; i < pop_size*prob_size; i++) {
        int rand_bit = (int)round(drand48()); // random 0 or 1
        pop[i] = rand_bit;
    }
}

/* One Max fitness function - count the number of ones in bitstring */
int assess_fitness(int* indiv, int prob_size) {
    /* TASK 1.1: IMPLEMENT FITNESS FUNCTION */
    int sum = 0;

    for (int i = 0; i < prob_size; i++) {
        //printf("indiv: %d, size: %d\n", indiv[i], prob_size);
        sum += indiv[i];
    }
    
    return sum;
}

// Compare function for qsort - sorts fitness in descending order
int compare(const void* a, const void* b) {
    const pop_fitness_map* a_cast  = a;
    const pop_fitness_map* b_cast  = b;
    if (a_cast->fitness < b_cast->fitness) {
        return 1;
    } else if (a_cast->fitness == b_cast->fitness) {
        return 0;
    } else {
        return -1;
    }
}

// Select the num_elites fittest indivs from the population
void select_elites(int* pop, int* elites, int pop_size, int prob_size, int num_elites) {
    pop_fitness_map pop_fitness_map_arr[pop_size];
    // populate struct to hold idx and fitness vals
    for (int i = 0; i < pop_size; i++) {
        int* indiv = &pop[i*prob_size];
        int fitness = assess_fitness(indiv, prob_size);
        pop_fitness_map_arr[i].pop_idx = i;
        pop_fitness_map_arr[i].fitness = fitness;
    }
    // do a sort on the struct to put fitness vals in descending order
    qsort(pop_fitness_map_arr, pop_size, sizeof(pop_fitness_map), compare);
    // extract elites via sorted idx
    for (int i = 0; i < num_elites; i++) {
        int pop_idx = pop_fitness_map_arr[i].pop_idx;
        int* elite = &pop[pop_idx*prob_size];
        memcpy(&elites[prob_size*i], elite, sizeof(int)*prob_size);
    }
}

/* Tournament selection (with tournament size of 2)
 * Algorithm 32 from 'Essentials of Metaheuristics'
 */
void select_parents(int* pop, int* parents, int pop_size, int prob_size) {
    for (int i = 0; i < 2; i++) { // 2 loops to select 2 parents
		// Pick an individual at random to be the best
		int rand_pop_index = (int)floor(drand48()*pop_size);
		int best_index = rand_pop_index;
        int* best_indiv = &pop[best_index*prob_size];
		int best_fitness = assess_fitness(best_indiv, prob_size);
        /* Using tournament size of 2, so select one more random indiv to
         * compete */
		for (int j = 1; j < 2; j++) {
			rand_pop_index = (int)floor(drand48()*pop_size);
			int* indiv = &pop[rand_pop_index*prob_size];
			int fitness = assess_fitness(indiv, prob_size);
			if (fitness > best_fitness) {
				best_indiv = indiv;
				best_index = rand_pop_index;
			}
		}
        // set parent
        memcpy(&parents[prob_size*i], best_indiv, sizeof(int)*prob_size);
	}
}

/* Uniform crossover
 * Algorithm 25 from 'Essentials of Metaheuristics'
 */
void crossover(int* parents, int* children, int prob_size) {
    // parents and children are arrays of length 2*prob_size
    /* TASK 1.2: SET CROSSOVER PROBABILITY */
    double pr_swap = 0.0;
    for (int i = 0; i < prob_size; i++) {
        double rand_double = drand48();
        if (rand_double <= pr_swap) {
            /* give child a value of parent b at this index and child b value
             of parent a */
            children[i] = parents[prob_size+i];
            children[prob_size+i] = parents[i];
        } else {
            // children get their parents' values at this index
            children[i] = parents[i];
            children[prob_size+i] = parents[prob_size+i];
        }
    }
}

/* Bit-flip mutation
 * Algorithm 22 from 'Essentials of Metaheuristics'
 */
void mutate(int* indiv, int* indiv_mut, int prob_size) {
    /* TASK 1.3: SET MUTATION PROBABILITY */
    double pr_flip = 0.0;
    for (int i = 0; i < prob_size; i++) {
        double rand_double = drand48();
        if (rand_double <= pr_flip) {
            // flip the bit
            if (indiv[i] == 0) {
                indiv_mut[i] = 1;
            } else {
                indiv_mut[i] = 0;
            }
        } else {
            // copy the bit
            indiv_mut[i] = indiv[i];
        }
    }
}

// Fill north and south send buffers with appropriate indivs from pop
void prepare_migrants(int* pop, int* north_send, int* south_send, int side_len,
        int pop_size, int prob_size) {
    /* TASK 2.5: CALCULATE THE SOUTH START INDEX */
    int south_start_index = 0;
    for(int i = 0; i < pop_size; i ++){
        if(pop[i] == 1){
            south_start_index = i;
            break;
        }
    }
    memcpy(north_send, &pop[0], sizeof(int)*side_len*prob_size);
    memcpy(south_send, &pop[south_start_index], sizeof(int)*side_len*prob_size);
}

/* Replace apprpriate indivs in pop with indivs from north and south recv
 * buffers
 */
void integrate_migrants(int* pop, int* north_recv, int* south_recv, int side_len,
        int pop_size, int prob_size) {
    /* TASK 2.5: CALCULATE THE SOUTH START INDEX */
    int south_start_index = 0;
    for(int i = 0; i < pop_size; i ++){
        if(pop[i] == 1){
            south_start_index = i;
            break;
        }
    }
    memcpy(&pop[0], north_recv, sizeof(int)*side_len*prob_size);
    memcpy(&pop[south_start_index], south_recv, sizeof(int)*side_len*prob_size);
}

void print_usage_info(void) {
    fprintf(stderr, "Usage: ga_island.exe pop_size prob_size seed\n");
    fprintf(stderr, "Mandatory arguments:\n");
    fprintf(stderr, "  pop_size: number of individuals in population "
            " of each rank (must be even square number >= 4)\n");
    fprintf(stderr, "  prob_size: size of One Max problem to solve\n");
    fprintf(stderr, "  seed: seed value for generating random numbers\n");
}

int main(int argc, char** argv) {
    // Setup MPI vars
    int my_rank, comm_size;
    const MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &comm_size);

    // Check and parse args
    if (argc != 4) {
        print_usage_info();
        exit(1);
    }

    // Enforce even square >= 4 local pop size
    const int pop_size = atoi(argv[1]);
    if (pop_size < 4) {
        fprintf(stderr, "Local population size must be even square number >= 4\n");
        exit(2);
    }
    double sr = sqrt(pop_size);
    int sr_int = (int)round(sr);
    if ((int)pow(sr_int, 2) != pop_size) {
        fprintf(stderr, "Local population size must be even square number >= 4\n");
        exit(3);
    }
    if (pop_size % 2 != 0) {
        fprintf(stderr, "Local population size must be even square number >= 4\n");
        exit(4);
    }
    const int side_len = sr_int; // length of square edge in local pop

    const int prob_size = atoi(argv[2]); // length of bitstrings

    // seed RNG dependent on rank
    const long int iseedlong_base = atoi(argv[3]);
    /* TASK 2.1: CHANGE THE SEED TO BE UNIQUE FOR EACH RANK */
    struct drand48_data* randBuffer;
    for(int i = 0; i < comm_size; i++){
        const long int iseedlong_rank = iseedlong_base + i;
        randBuffer = malloc(sizeof(struct drand48_data));
        srand48_r(iseedlong_rank, randBuffer);
    }
    //srand48(iseedlong_rank);

    root_print(my_rank, "Using local population size of %d\n", pop_size);
    root_print(my_rank, "Global population size is %d over %d ranks\n",
            pop_size*comm_size, comm_size);
    root_print(my_rank, "Solving One Max problem of size %d\n", prob_size);
    root_print(my_rank, "Using base RNG seed of %ld\n", iseedlong_base);

    // init local population
    int* current_pop = (int*)malloc(sizeof(int)*pop_size*prob_size);
    init_population(current_pop, pop_size, prob_size);

    // Main generational loop
    int num_gens = 1;
    bool found_ideal = false;
    // set num elites so it's always even but proportional to pop_size, roughly
    // 20% of pop size
    const int num_elites = (int)(round(pop_size*0.1)*2);
    root_print(my_rank, "Number of elites to preserve between each generation: %d\n", num_elites);
    while (true) {
        root_print(my_rank, "\nGeneration %d\n", num_gens);
        // Assess fitness of each individual in the population
        int best_fitness = 0;
        for (int i = 0; i < pop_size; i++) {
            int* indiv = &current_pop[i*prob_size];
            int fitness = assess_fitness(indiv, prob_size);
            // update best fitness for this generation if necessary
            if (fitness > best_fitness) {
                best_fitness = fitness;
            }
        }

        // Determine if an ideal individual exists across the ranks
        /* Make buffer to receive best fitness vals (only actually used
         * on the root rank)
         */
        int* best_fitness_vals = (int*)malloc(sizeof(int)*comm_size);
        // Receive best fitness vals from all ranks on root rank
        MPI_Gather(&best_fitness, 1, MPI_INT, best_fitness_vals, 1,
                MPI_INT, root, comm);
        if (my_rank == root) {
            int best_fitness_this_gen = 0;
            int best_rank = 0;
            for (int i = 0; i < comm_size; i++) {
                int fitness = best_fitness_vals[i];
                if (fitness > best_fitness_this_gen) {
                    best_fitness_this_gen = fitness;
                    best_rank = i;
                }
                if (fitness == prob_size) {
                    found_ideal = true;
                    printf("Ideal solution found on rank %d\n", i);
                    printf("Sending termination signal to all ranks then "
                            "exiting...\n");
                    break;
                }
            }
            printf("Best fitness = %d, on rank %d\n", best_fitness_this_gen,
                best_rank);
        }
        free(best_fitness_vals);

        // Check if it is time to terminate
        /* TASK 2.6: ADD TERMINATION CONDITION FOR MAIN WHILE LOOP */
        if(found_ideal){
            MPI_Bcast(&comm_size, comm_size, MPI_INT, root, comm);
            break;
        }

        // Migration between islands
        // Construct buffers for sending/receiving
        int* north_recv = (int*)malloc(sizeof(int)*side_len*prob_size);
        int* south_recv = (int*)malloc(sizeof(int)*side_len*prob_size);
        int* north_send = (int*)malloc(sizeof(int)*side_len*prob_size);
        int* south_send = (int*)malloc(sizeof(int)*side_len*prob_size);
        MPI_Request recv_reqs[2];
        MPI_Status recv_statuses[2];
        MPI_Request send_reqs[2];
        MPI_Status send_statuses[2];
        /* Communication pattern is:
         * recv N
         * recv S
         * send N
         * send S
         * wait on data received from N and S
         * select elites and do genetic operations
         * wait on data sent to N and S
         * update population
         */
        /* TASK 2.2: DETERMINE rank_north AND rank_south VALUES */
        int rank_north = my_rank - 1;
        int rank_south = my_rank + 1;

        /* TASK 2.3: RECEIVE MIGRANTS */
        MPI_Request request;
        MPI_Irecv(&rank_north, 1, MPI_INT, rank_north, rank_north, MPI_COMM_WORLD, &request); // north recv
        MPI_Irecv(&rank_south, 1, MPI_INT, rank_north, rank_south, MPI_COMM_WORLD, &request); // south recv
        
        // fill north and south send buffers with migrant indivs
        prepare_migrants(current_pop, north_send, south_send, side_len,
                pop_size, prob_size);

        /* TASK 2.4: SEND MIGRANTS */
        MPI_Isend(&rank_north, 1, MPI_INT, rank_north, rank_north, MPI_COMM_WORLD, &request); // north send
        MPI_Isend(&rank_south, 1, MPI_INT, rank_south, rank_south, MPI_COMM_WORLD, &request); // south send

        // Wait on both receives
        MPI_Waitall(2, recv_reqs, recv_statuses);

        // Integrate received migrants into population
        integrate_migrants(current_pop, north_recv, south_recv, side_len,
                pop_size, prob_size);
        free(north_recv);
        free(south_recv);
        free(north_send);
        free(south_send);

        // init new pop: set Q
        int* new_pop = (int*)malloc(sizeof(int)*pop_size*prob_size);
        // Select elites from pop
        int elites[num_elites*prob_size];
        select_elites(current_pop, elites, pop_size, prob_size, num_elites);

        // Create new local population via genetic operations
        for (int i = 0; i < (pop_size - num_elites); i+=2) {
            int parents[2*prob_size];
            select_parents(current_pop, parents, pop_size, prob_size);
            int children[2*prob_size];
			crossover(parents, children, prob_size);
            // mutate both children
            int* child_one = &children[0];
            int* child_two = &children[prob_size];
            int children_mut[2*prob_size];
            int* child_one_mut = &children_mut[0];
            int* child_two_mut = &children_mut[prob_size];
            mutate(child_one, child_one_mut, prob_size);
            mutate(child_two, child_two_mut, prob_size);
            // add children into new population
            memcpy(&new_pop[i*prob_size], child_one_mut, sizeof(int)*prob_size);
            memcpy(&new_pop[(i+1)*prob_size], child_two_mut,
                    sizeof(int)*prob_size);
        }
        // Wait on both sends
        MPI_Waitall(2, send_reqs, send_statuses);

        // add elites to new population (at end) - this ensures they don't
        // accidentally get removed in above for loop
        memcpy(&new_pop[(pop_size - num_elites)*prob_size], elites, sizeof(int)*num_elites*prob_size);
        // Update population: P <- Q
        memcpy(current_pop, new_pop, sizeof(int)*pop_size*prob_size);
		num_gens++;
        free(new_pop);
    }
    free(current_pop);
    free(randBuffer); 
    root_print(my_rank, "Found ideal solution after %d generations."
        " Exiting...\n", num_gens);
    MPI_Finalize();
    return 0;
}
