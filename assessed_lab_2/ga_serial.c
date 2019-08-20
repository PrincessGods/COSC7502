/*
 * ga_serial.c by Jordan T. Bishop 2019
 * Developed for COSC3500/7502 at the University of Queensland
 *
 * Compile with:
 * gcc -std=gnu99 ga_serial.c -lm -o ga_serial.exe -O3 -Wall -Wextra
 *
 * Run with ga_serial.sh Slurm submission script.
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <stdarg.h>

// struct for use with qsort comparator
typedef struct {
    int pop_idx;
    int fitness;
} pop_fitness_map;

bool verbose = false;

void verbose_print(char* format, ...) {
    va_list args;
    if (!verbose) {
        return;
    }
    va_start(args, format);
    vprintf(format, args);
    va_end(args);
}

// Print single indiv bitstring to stdout
void print_indiv(int* indiv, int prob_size) {
   for (int i = 0; i < prob_size; i++) {
        verbose_print("%d", indiv[i]);
   }
}

// Print genetic operations to stdout
void print_gen_ops(int* parents, int* children, int* children_mut,
        int prob_size) {
    print_indiv(&parents[0], prob_size);
    verbose_print(" X ");
    print_indiv(&parents[prob_size], prob_size);
    verbose_print(" -> ");
    print_indiv(&children[0], prob_size);
    verbose_print(", ");
    print_indiv(&children[prob_size], prob_size);
    verbose_print(" ~ ");
    print_indiv(&children_mut[0], prob_size);
    verbose_print(", ");
    print_indiv(&children_mut[prob_size], prob_size);
    verbose_print("\n");
}

void print_population(int* pop, int pop_size, int prob_size) {
    for (int i = 0; i < pop_size; i++) {
        verbose_print("Individual %d: ", i+1);
        for (int j = 0; j < prob_size; j++) {
            verbose_print("%d", pop[i*prob_size + j]);
        }
        verbose_print("\n");
    }
}

void print_elites(int* elites, int num_elites, int prob_size) {
    for (int i = 0; i < num_elites; i++) {
        verbose_print("Elite %d: ", i+1);
        for (int j = 0; j < prob_size; j++) {
            verbose_print("%d", elites[i*prob_size + j]);
        }
        verbose_print("\n");
    }
}

/* Fill the pop array with random zeros/ones
 * Algorithm 21 from 'Essentials of Metaheuristics'
 */
void init_population(int* pop, int pop_size, int prob_size) {
    printf("Initialising population\n");
    for (int i = 0; i < pop_size*prob_size; i++) {
        int rand_bit = (int)round(drand48()); // random 0 or 1
        pop[i] = rand_bit;
    }
}

int assess_fitness(int* indiv, int prob_size) {
    /* TASK 1.1: IMPLEMENT FITNESS FUNCTION */
    int count = 0;
    int sum = 0;
    while(count <= prob_size){
        if(indiv == 1){
            sum ++;
        }

        count ++
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
    if (num_elites > 0) {
        verbose_print("Selecting elites\n");
    }
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

void print_usage_info(void) {
    fprintf(stderr, "Usage: ga_serial.exe pop_size prob_size seed [-v]\n");
    fprintf(stderr, "Mandatory arguments:\n");
    fprintf(stderr, "  pop_size: number of individuals in population "
            "(must be even)\n");
    fprintf(stderr, "  prob_size: size of One Max problem to solve\n");
    fprintf(stderr, "  seed: seed value for generating random numbers\n");
    fprintf(stderr, "Optional arguments:\n");
    fprintf(stderr, "  -v: enable verbose output "
            "(prints all generational info)\n");
}

int main(int argc, char** argv) {
    // Check and parse args
    if (argc != 4) {
        if (argc == 5) {
            if (strcmp(argv[4], "-v") != 0) {
                print_usage_info();
                exit(1);
            } else {
                verbose = true;
            }
        } else {
            print_usage_info();
            exit(2);
        }
    }
    const int pop_size = atoi(argv[1]);
    // Enforce minimum pop size that is even
    if (pop_size < 4 || pop_size % 2 != 0) {
        printf("Population size must be >=4 and even! "
                "Specified population size was %d\n", pop_size);
        exit(3);
    }
    const int prob_size = atoi(argv[2]); // length of bitstrings
    int* current_pop = (int*)malloc(sizeof(int)*pop_size*prob_size);
    // seed RNG
    const long int iseedlong = atoi(argv[3]);
    srand48(iseedlong);

    printf("Using population size of %d\n", pop_size);
    printf("Solving One Max problem of size %d\n", prob_size);
    printf("Using seed of %ld\n", iseedlong);
    verbose_print("Verbose output enabled\n");

    // init local population
    init_population(current_pop, pop_size, prob_size);

    // Main generational loop
    int num_gens = 1;
    // set num elites so it's always even but proportional to pop_size, roughly
    // 20% of pop size
    const int num_elites = (int)(round(pop_size*0.1)*2);
    printf("Number of elites to preserve between each generation: %d\n", num_elites);
    while (true) {
        printf("\nGeneration %d\n", num_gens);
        print_population(current_pop, pop_size, prob_size);
        // Assess fitness of each individual in the population
        int best_fitness = 0;
        for (int i = 0; i < pop_size; i++) {
            int* indiv = &current_pop[i*prob_size];
            int fitness = assess_fitness(indiv, prob_size);
            verbose_print("Fitness for individual %d = %d\n", i+1, fitness);
            // update best fitness for this generation if necessary
            if (fitness > best_fitness) {
                best_fitness = fitness;
            }
        }

        printf("Best fitness = %d\n", best_fitness);
        // Determine if an ideal individual has been found
        if (best_fitness == prob_size) {
            break;
        }

        // init new pop: set Q
        int* new_pop = (int*)malloc(sizeof(int)*pop_size*prob_size);
        // Select elites from pop
        int elites[num_elites*prob_size];
        select_elites(current_pop, elites, pop_size, prob_size, num_elites);
        print_elites(elites, num_elites, prob_size);

        // Create new local population via genetic operations
        verbose_print("Applying genetic operations\n");
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
            print_gen_ops(parents, children, children_mut, prob_size);
            // add children into new population
            memcpy(&new_pop[i*prob_size], child_one_mut, sizeof(int)*prob_size);
            memcpy(&new_pop[(i+1)*prob_size], child_two_mut,
                    sizeof(int)*prob_size);
        }
        // add elites to new population (at end) - this ensures they don't
        // accidentally get removed in above for loop
        memcpy(&new_pop[(pop_size - num_elites)*prob_size], elites, sizeof(int)*num_elites*prob_size);
        // Update population: P <- Q
        memcpy(current_pop, new_pop, sizeof(int)*pop_size*prob_size);
		num_gens++;
        free(new_pop);
    }
    free(current_pop);
    printf("Found ideal solution after %d generations. Exiting...\n", num_gens);
    return 0;
}
