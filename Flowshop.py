#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import numpy as np
import copy

'''==================== main code ==============================='''
'''----- generate initial population -----'''
def geneAlgo(pt, population_size, crossover_rate, mutation_rate, num_mutation_jobs, num_iteration, num_job, num_m):
    Tbest=999999999999999
    best_list,best_obj=[],[]
    population_list=[]
    for i in range(population_size):
        random_num=list(np.random.permutation(num_job)) # generate a random permutation of 0 to num_job-1
        population_list.append(random_num) # add to the population_list

    for n in range(num_iteration):
        Tbest_now=99999999999           
        '''-------- crossover --------'''
        parent_list=copy.deepcopy(population_list)
        offspring_list=copy.deepcopy(population_list)
        S=list(np.random.permutation(population_size)) # generate a random sequence to select the parent chromosome to crossover

        for m in range(int(population_size/2)):
            crossover_prob=np.random.rand()
            if crossover_rate>=crossover_prob:
                parent_1= population_list[S[2*m]][:]
                parent_2= population_list[S[2*m+1]][:]
                child_1=['na' for i in range(num_job)]
                child_2=['na' for i in range(num_job)]
                fix_num=round(num_job/2)
                g_fix=list(np.random.choice(num_job, fix_num, replace=False))

                for g in range(fix_num):
                    child_1[g_fix[g]]=parent_2[g_fix[g]]
                    child_2[g_fix[g]]=parent_1[g_fix[g]]
                c1=[parent_1[i] for i in range(num_job) if parent_1[i] not in child_1]
                c2=[parent_2[i] for i in range(num_job) if parent_2[i] not in child_2]

                for i in range(num_job-fix_num):
                    child_1[child_1.index('na')]=c1[i]
                    child_2[child_2.index('na')]=c2[i]
                offspring_list[S[2*m]]=child_1[:]
                offspring_list[S[2*m+1]]=child_2[:]

        '''--------mutation--------'''   
        for m in range(len(offspring_list)):
            mutation_prob=np.random.rand()
            if mutation_rate >= mutation_prob:
                m_chg=list(np.random.choice(num_job, num_mutation_jobs, replace=False)) # chooses the position to mutation
                t_value_last=offspring_list[m][m_chg[0]] # save the value which is on the first mutation position
                for i in range(num_mutation_jobs-1):
                    offspring_list[m][m_chg[i]]=offspring_list[m][m_chg[i+1]] # displacement

                offspring_list[m][m_chg[num_mutation_jobs-1]]=t_value_last # move the value of the first mutation position to the last mutation position


        '''--------fitness value(calculate idle time)-------------'''
        total_chromosome=copy.deepcopy(parent_list)+copy.deepcopy(offspring_list) # parent and offspring chromosomes combination
        chrom_fitness,chrom_fit=[],[]
        total_fitness=0
        s,d,D={},{},{}
        v=[0 for c in range(population_size*2)]

        for c in range(population_size*2):
            for i in range(num_m):
                s[c,i]=pt[total_chromosome[c][0]][i]
                d[c,i]=v[c]
                D[c,i,total_chromosome[c][0]]=v[c]
                v[c]=v[c]+pt[total_chromosome[c][0]][i]

            for j in range(num_job):
                D[c,0,j]=0

            for j in range(1,num_job):
                for i in range(0,num_m-1):
                    s[c,i]=s[c,i]+pt[total_chromosome[c][j]][i]
                    D[c,i+1,j]=max(0,s[c,i]+d[c,i]-s[c,i+1]-d[c,i+1])
                    d[c,i+1]=d[c,i+1]+D[c,i+1,j]

                s[c,num_m-1]=s[c,num_m-1]+pt[total_chromosome[c][j]][i+1]

            v[c]=0
            for i in range(num_m):
                v[c]=v[c]+d[c,i]

            chrom_fitness.append(1/v[c])
            chrom_fit.append(v[c])
            total_fitness=total_fitness+chrom_fitness[c]

        '''----------selection----------'''
        pk,qk=[],[]

        for i in range(population_size*2):
            pk.append(chrom_fitness[i]/total_fitness)
        for i in range(population_size*2):
            cumulative=0
            for j in range(0,i+1):
                cumulative=cumulative+pk[j]
            qk.append(cumulative)

        selection_rand=[np.random.rand() for i in range(population_size)]

        for i in range(population_size):
            if selection_rand[i]<=qk[0]:
                population_list[i]=copy.deepcopy(total_chromosome[0])
            else:
                for j in range(0,population_size*2-1):
                    if selection_rand[i]>qk[j] and selection_rand[i]<=qk[j+1]:
                        population_list[i]=copy.deepcopy(total_chromosome[j+1])
                        break

        '''----------comparison----------'''
        for i in range(population_size*2):
            if chrom_fit[i]<Tbest_now:
                Tbest_now=chrom_fit[i]
                sequence_now=copy.deepcopy(total_chromosome[i])

        if Tbest_now<=Tbest:
            Tbest=Tbest_now
            sequence_best=copy.deepcopy(sequence_now)
            
    return sequence_best, Tbest

