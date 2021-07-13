
/*
	JADE.java
	Isaac Triguero Velazquez.
	
	Created by Isaac Triguero Velazquez  23-7-2009
	Copyright (c) 2008 __MyCompanyName__. All rights reserved.

*/

package keel.Algorithms.Instance_Generation.LSHADE;

import keel.Algorithms.Instance_Generation.Basic.PrototypeSet;
import keel.Algorithms.Instance_Generation.Basic.PublicVar;
import keel.Algorithms.Instance_Generation.Basic.PrototypeGenerator;
import keel.Algorithms.Instance_Generation.Basic.Prototype;
import keel.Algorithms.Instance_Generation.Basic.PrototypeGenerationAlgorithm;
import keel.Algorithms.Instance_Generation.Chen.ChenGenerator;
import keel.Algorithms.Instance_Generation.HYB.HYBGenerator;
import keel.Algorithms.Instance_Generation.*;
import java.util.*;

import keel.Algorithms.Instance_Generation.utilities.*;
import keel.Algorithms.Instance_Generation.utilities.KNN.*;

import org.core.*;

import org.core.*;


import java.util.StringTokenizer;




/**
 * @param k Number of neighbors
 * @param Population Size.
 * @param ParticleSize.
 * @param Scaling Factor.
 * @param Crossover rate.
 * @param Strategy (1-5).
 * @param MaxIter
 * @author Isaac Triguero
 * @version 1.0
 */
public class LSHADEGenerator extends PrototypeGenerator {

  /*Own parameters of the algorithm*/
  
  // We need the variable K to use with k-NN rule
  private int k;
 
  private int PopulationSize; 
  private int ParticleSize;
  private int MaxIter; 
  private double ScalingFactor;
  private double CrossOverRate[];
  private int Strategy;
  private String CrossoverType; // Binomial, Exponential, Arithmetic
  //private boolean Archive; // JAde with or without Archive.
  private double arc_rate; // to select the number of pbest
  private double c;
  protected int numberOfClass; 
  protected int numberOfbetters; // numero de mejores atener en cuenta
  protected int numberOfPrototypes;  // Particle size is the percentage
  /** Parameters of the initial reduction process. */
  private String[] paramsOfInitialReducction = null;
  private int MaxEval;
  
  /**
   * Build a new JADEGenerator Algorithm
   * @param t Original prototype set to be reduced.
   * @param perc Reduction percentage of the prototype set.
   */
  
  public LSHADEGenerator(PrototypeSet _trainingDataSet, int neigbors,int poblacion, int perc, int iteraciones, double F, double CR, int strg)
  {
      super(_trainingDataSet);
      algorithmName="JADE";
      
      this.k = neigbors;
      this.PopulationSize = poblacion;
      this.ParticleSize = perc;
      this.MaxIter = iteraciones;
      this.numberOfPrototypes = getSetSizeFromPercentage(perc);
      
      this.ScalingFactor = F;
      this.MaxEval = (int)(PublicVar.EvalWeight*trainingDataSet.get(0).numberOfInputs()*trainingDataSet.size());
      this.Strategy = strg;
      
  }
  


  /**
   * Build a new JADEGenerator Algorithm
   * @param t Original prototype set to be reduced.
   * @param params Parameters of the algorithm (only % of reduced set).
   */
  public LSHADEGenerator(PrototypeSet t, Parameters parameters)
  {
      super(t, parameters);
      algorithmName="JADE";
      this.k =  parameters.getNextAsInt();
      this.PopulationSize =  parameters.getNextAsInt();
      this.ParticleSize =  parameters.getNextAsInt();
      this.MaxIter =  parameters.getNextAsInt();
      this.arc_rate = parameters.getNextAsDouble();
      this.arc_rate =1.4; //Copy optimal value from the author. 
      this.c =parameters.getNextAsDouble();
       
      this.numberOfPrototypes = getSetSizeFromPercentage(ParticleSize);
      this.numberOfbetters= (int) (this.arc_rate*PopulationSize);
      if( numberOfbetters <1) numberOfbetters  = 1;
      this.MaxEval = (int)(PublicVar.EvalWeight*trainingDataSet.get(0).numberOfInputs()*trainingDataSet.size());
      System.out.println("Numero de p-best = "+ this.numberOfbetters);
      this.numberOfClass = trainingDataSet.getPosibleValuesOfOutput().size();
      System.out.print("\nIsaac dice:  " + k + " Swar= "+PopulationSize+ " Particle=  "+ ParticleSize + " Maxiter= "+ MaxIter+" CR=  "+this.CrossOverRate+ "\n");
      //numberOfPrototypes = getSetSizeFromPercentage(parameters.getNextAsDouble());
  }
  
  //mejoresParticulas: best particular
  //Psxhll: Return the index array  meaning the descending order of fitness. 
  public int[] mejoresParticulas(PrototypeSet population[],double fitness[]){
	//int number = this.numberOfbetters;
	//  int index[] = new int[number];
	  int index[] = new int[population.length];
	  int ind= 0;
	  double mejor = Double.MIN_VALUE;
	  double acc;
	  
	  for(int i=0; i< population.length; i++){
		  index[i] = -1; //Psxhll: Init value -1 to avoid index 0
		  acc =fitness[i]; //accuracy(population[i],trainingDataSet); 
		  
		  if(acc > mejor )
		  {
			ind = i;
			mejor = acc;
		  }	  
	  }
	  index[0] = ind;
	  
	  for (int j=1; j<population.length; j++){
		  mejor = Double.MIN_VALUE;
		  for(int i=0; i< population.length; i++){
			  acc = fitness[i];//accuracy(population[i],trainingDataSet); 
			  
			  //if(acc > mejor && acc < accuracy(population[index[j-1]],trainingDataSet))
			  if((acc > mejor && acc <= fitness[index[j-1]])) { 
					 // (acc is not in index) //Psxhll How about two same fitness values?
				if (!isExist(index, i)) { 
					if (i == 14) {
						int a = 1;
					}
					ind = i;
					mejor = acc;
				}
			  }	  
		  }
		  index[j] =  ind;
	  
	  }
	  return index;
  }
  
  private boolean isExist(int[]index, int pos) {
	  boolean result = false;
	  for (int i = 0; i< index.length; i++) {
		  if (index[i] == pos) {
			  result = true;
			  break;
		  }
	  }
	  return result;
  }
  public PrototypeSet mutant(PrototypeSet population[], double fitness[], int actual,  PrototypeSet Archivo[],int utilArchivo){
	  
	  
	  PrototypeSet mutant = new PrototypeSet(population.length);
	  PrototypeSet r1,r2,xbest, resta, producto, resta2, producto2, result;
	  
		   
	  // r1 different to actual	      
	  int ran;
	   do{
		   ran =  RandomGenerator.Randint(0, population.length);   
	   }while(ran == actual);
	   
	   r1 = population[ran];
	   
	  int number;
	  
	  do{ 
		  number = RandomGenerator.Randint(0, population.length+ utilArchivo);
	  }while (number==ran || number == actual );
	  
	   if(number < population.length){
		 r2 = population[number];  
	   }else
		 r2 = Archivo[number-population.length];
	   
	   // Tengo que sacar los 100p % mejores de la poblaci�n actual.
	  // System.out.println("Numero de p-best = "+ num_mejores);
	  
	   
	   int indices[] = new int [this.numberOfbetters];
	   
	   indices = mejoresParticulas(population,fitness);
	   number = RandomGenerator.Randint(0, indices.length);
	   
	   xbest = population[indices[number]];

	   /*
	   for(int i=0; i< population.length; i++){
		   System.out.println(accuracy(population[i],trainingDataSet));
	   }
	   for( int i=0; i< num_mejores ; i++){
		   System.out.println(indices[i]);
	   }
	   */
			switch(this.Strategy){
		   	   case 1: 
		   		   resta = xbest.restar(population[actual]);
		   		   resta2 = r1.restar(r2);
		   		   
		   		   producto = resta.mulEscalar(this.ScalingFactor);
		   		   producto2 = resta2.mulEscalar(this.ScalingFactor);
		   		   
		   		   result = producto.sumar(producto2);
		   		   mutant = population[actual].sumar(result);
		   	    break;
			   		   	
		   }   
	   

	  // System.out.println("********Mutante**********");
	 // mutant.print();
	   
     mutant.applyThresholds();
	
	  return mutant;
  }
  
	  /*
	  Return random value from Cauchy distribution with mean "mu" and variance "gamma"
	  http://www.sat.t.u-tokyo.ac.jp/~omi/random_variables_generation.html#Cauchy
	*/
  private double cauchy_g(double mu, double gamma) {
	return mu + gamma * Math.tan(Math.PI*(Math.random() - 0.5));
	}
	
	/*
	  Return random value from normal distribution with mean "mu" and variance "gamma"
	  http://www.sat.t.u-tokyo.ac.jp/~omi/random_variables_generation.html#Gauss
	*/
	private double gauss(double mu, double sigma){
	return mu + sigma * Math.sqrt(-2.0 * Math.log(Math.random())) * Math.sin(2.0 * Math.PI * Math.random());
	}

  /**
   * Generate a reduced prototype set by the JADEGenerator method.
   * @return Reduced set by JADEGenerator's method.
   */
  
  
  public PrototypeSet reduceSet()
  {
	  System.out.print("\nThe algorithm  JADE is starting...\n Computing...\n");
	  
	  System.out.println("Number of prototypes, result set = "+numberOfPrototypes+ "\n");
	  
	  if(numberOfPrototypes < trainingDataSet.getPosibleValuesOfOutput().size()){
		  System.out.println("Number of prototypes less than the number of clases");
		  numberOfPrototypes = trainingDataSet.getPosibleValuesOfOutput().size();
	  }
	  System.out.println("Reduction %, result set = "+((trainingDataSet.size()-numberOfPrototypes)*100)/trainingDataSet.size()+ "\n");

	  
	//  System.out.println("training Size->" +trainingDataSet.size());
	  
	  //Algorithm
	  // First, we create the population, with PopulationSize.
	  // like a prototypeSet's vector.  
	  PrototypeSet nominalPopulation;
	  PrototypeSet population [] = new PrototypeSet [PopulationSize];
	  PrototypeSet mutation[] = new PrototypeSet[PopulationSize];
	  PrototypeSet crossover[] = new PrototypeSet[PopulationSize];
	  double[] dif_fitness = new double[PopulationSize];
	//Psxhll: Memories of sf and cr have 5 elements, help generata new cr and sf at a new generation. They are computed at the end of a for-loop based on
      //fitness_diff[] and all sucessfull values of cr and sf. The memory size is limited, predefied by users.
	  int memory_size = 5;
		double[] memory_sf = new double[memory_size];
		double[] memory_cr = new double[memory_size];
	
		for (int i = 0; i < memory_size; i++) {
		    memory_sf[i] = 0.5;
		    memory_cr[i] = 0.5;
		}

	//memory index counter
	
	  int memory_pos =0;
	  double fitness[] = new double[PopulationSize];
	  double arc_rate = 1.4;
	  double meanCR = 0.5;
	  double meanF = 0.5;
	  int num_arc_inds =0;
	  //Psxhll: archivo store all parents that gennerate better offspring.  
	  PrototypeSet Archivo[] = new PrototypeSet[ (int)(arc_rate*this.PopulationSize)];
	  int utilArchivo = 0;
	  
	  // we save the differents successful F and CR.
	  double SF[] = new double[this.PopulationSize];
	  double SCR[] = new double[this.PopulationSize];
	  
	  this.CrossOverRate= new double[this.PopulationSize];
	  double F[] = new double[this.PopulationSize];
	  

	  // First Stage, Initialization.
	  population[0]=selecRandomSet(numberOfPrototypes,true).clone() ;
	  
	  // Aseguro que al menos hay un representante de cada clase.
	  PrototypeSet clases[] = new PrototypeSet [this.numberOfClass];
	  for(int i=0; i< this.numberOfClass; i++){
		  clases[i] = new PrototypeSet(trainingDataSet.getFromClass(i));
	  }
	
	  for(int i=0; i< population[0].size(); i++){
		  for(int j=0; j< this.numberOfClass; j++){
			  if(population[0].getFromClass(j).size() ==0 && clases[j].size()!=0){
				  
				  population[0].add(clases[j].getRandom());
			  }
		  }
	  }
	  
	  //population[0].print();
	  
	  fitness[0] = accuracy(population[0],trainingDataSet); 
	  
	 // population[0].print();
	  for(int i=1; i< PopulationSize; i++){
		  population[i] = new PrototypeSet();
		  for(int j=0; j< population[0].size(); j++){
			  population[i].add(trainingDataSet.getFromClass(population[0].get(j).getOutput(0)).getRandom());
		  }
		  fitness[i] = accuracy(population[i],trainingDataSet);   // DE fitness, no hace falta formatear porque son aleatorios!
	  }
	  
	  
	  //We select the best initial  particle
	 double bestFitness=fitness[0];
	  int bestFitnessIndex=0;
	  for(int i=1; i< PopulationSize;i++){
		  if(fitness[i]>bestFitness){
			  bestFitness = fitness[i];
			  bestFitnessIndex=i;
		  }
		  
	  }
	   for(int j=0;j<PopulationSize;j++){
         //Now, I establish the index of each prototype.
		  for(int i=0; i<population[j].size(); ++i)
			  population[j].get(i).setIndex(i);
	   }

	   System.out.println("Initial Fitness=  "+bestFitness);
	   this.Strategy = 1;
	   int ev = 0;
	   int LinearPopSize = PopulationSize;
	   double p_best_rate = 0.11;
	   
	   int p_num = (int)Math.round(PopulationSize *  p_best_rate);
	   int arc_size = (int)Math.round(PopulationSize * arc_rate);
	   while (ev < MaxEval)// Main loop
	   {  
			  int utilCR = 0;
			  // If we are going to use exponential, I calculate the index of possible selecting  Mutation.
			  for(int i=0; i<LinearPopSize; i++){
				  
				  
				  int rand_mem_index = (int)(Math.random()*memory_size);
				  double mu_sf = memory_sf[rand_mem_index];
				  double mu_cr = memory_cr[rand_mem_index];

					//generate CR_i and repair its value
					if (mu_cr == -1) {
						CrossOverRate[i] = 0;
					}
					else {
						CrossOverRate[i] = gauss(mu_cr, 0.1);
					    if (CrossOverRate[i] > 1) CrossOverRate[i] = 1;
					    else if (CrossOverRate[i] < 0) CrossOverRate[i] = 0;	
					}
			      
					//generate F_i and repair its value
					do {	
					    F[i] = cauchy_g(mu_sf, 0.1);
					} while (F[i] <= 0);
					if (F[i] > 1) 
						F[i] = 1;
				  
				 				  
				  this.ScalingFactor = F[i];
				  // Randomly choose xbestp, as one of the 100p% best Vector
				  
				   //Second:  Mutation Operation.
				   // I want to generate a PrototypeSet Mutation for each item of the population.
				   
				   mutation[i] = new PrototypeSet(population[i].clone());
			   
				   // Pasamos la poblaci�n, y la mejor actual, por si se usa /best/
				   mutation[i] = mutant(population,fitness, i, Archivo, utilArchivo).clone();

				   // Third: Crossver Operation.
				   // Now, we decide if the mutation will be the trial vector.

				   crossover[i] = new PrototypeSet(population[i].clone());
				   
				  for(int j=0; j< population[i].size(); j++){ // For each part of the solution
					   double randNumber = RandomGenerator.Randdouble(0, 1);
					   if(randNumber<this.CrossOverRate[i]){
						   crossover[i].set(j, mutation[i].get(j)); // Overwrite.
					   }
				   }
				   
				   // Fourth: Selection Operation.
				   // Decide if the trial vector is better than initial population.
				   //Crossover has the trialVector, we check its fitness.
				   
				  // crossover[i].applyThresholds();
				   nominalPopulation = new PrototypeSet();
			       nominalPopulation.formatear(population[i]);
			       fitness[i] = accuracy(nominalPopulation,trainingDataSet);
			       ev++;	  
			       nominalPopulation = new PrototypeSet();
			       nominalPopulation.formatear(crossover[i]);
			       
				  double trialVector = accuracy(nominalPopulation,trainingDataSet);
				  //System.out.println("Trial Vector fitness = "+ trialVector);
				  //System.out.println("fitness de la particula = "+ fitness[i]);
				  
				  if(trialVector > fitness[i]){
					  //System.out.println("Selecting");
					 
					  Archivo[utilArchivo%arc_size] = new PrototypeSet(population[i].clone());
					  utilArchivo++;
					  //successful parameters are preserved in S_F and S_CR
					  SCR[utilCR%PopulationSize] = this.CrossOverRate[i];
					  SF[utilCR%PopulationSize] = F[i];
					  dif_fitness[utilCR] = Math.abs(fitness[i] - trialVector);
					  utilCR++; //Psxhll: index of samples that we are collecting in the archive
					  population[i] = new PrototypeSet(crossover[i].clone());
					  fitness[i]  = trialVector;
					  utilArchivo = utilArchivo%arc_size;
				  }
				  if(fitness[i]>bestFitness){
					  bestFitness = fitness[i];
					  bestFitnessIndex=i;
					  System.out.println("Iter="+ ev +" Acc= "+ bestFitness);
				  }
			   }
			  // Now we remove solutions from A.
			  
			  if(utilArchivo > LinearPopSize){
				  utilArchivo = LinearPopSize;				  
			  }

			  //Psxhll: Update value of SCR and SF
			  
			  if (utilCR > 0) {
					double temp_sum_sf1 = 0;
					double temp_sum_sf2 = 0;
					double temp_sum_cr1 = 0;
					double temp_sum_cr2 = 0;
					double temp_sum = 0;
					double temp_weight = 0;

					for (int i = 0; i < utilCR; i++) 
						temp_sum += dif_fitness[i];
			      
					//weighted lehmer mean
					for (int i = 0; i < utilCR; i++) {
					    temp_weight = dif_fitness[i] / temp_sum;

					    temp_sum_sf1 += temp_weight * SF[i] * SF[i];
					    temp_sum_sf2 += temp_weight * SF[i];

					    temp_sum_cr1 += temp_weight * SCR[i] * SCR[i];
					    temp_sum_cr2 += temp_weight * SCR[i];
					}

					memory_sf[memory_pos] = temp_sum_sf1 / temp_sum_sf2;

					if (temp_sum_cr2 == 0 || memory_cr[memory_pos] == -1)
						memory_cr[memory_pos] = -1;
					else memory_cr[memory_pos] = temp_sum_cr1 / temp_sum_cr2;
					//increment the counter
					memory_pos++;
					if (memory_pos >= memory_size) memory_pos = 0;

			  }
			  

			 // System.out.println("Acc= "+ bestFitness);
			  //Psxhll: Reduce pop size
			  int num_redu_inds;
			  int min_pop_size = 4; //Psxhll: To be consistent with the param the author said in shade. 
			  int next_pop_size = (int)Math.round((((min_pop_size - PopulationSize) / (double)MaxIter) * ev) + PopulationSize);
			  //next_pop_size = (int)Math.round((((min_pop_size - initial_pop_size) / (double)max_func_evals) * nfes) + initial_pop_size);
			    if (LinearPopSize > next_pop_size) {
					num_redu_inds = LinearPopSize - next_pop_size;
					if ((LinearPopSize - num_redu_inds) <  min_pop_size) 
						num_redu_inds = LinearPopSize - min_pop_size;
	
					
					//we need to remove num_redu_inds instances from the population and fitness array. They are the worst in the pop
					int worst_ind;
					for (int k = 0; k < num_redu_inds; k++) {
					    worst_ind = 0;
					    for (int j = 1; j < LinearPopSize; j++) {
					    	if (fitness[j] < fitness[worst_ind]) 
					    		worst_ind = j;
					    }
	
					    //	    System.out.print("\ndeleted individual = " + worst_ind + "\n");
					    population = resizeMatrix(population, worst_ind);
					    //Archivo = resizeMatrix(Archivo, worst_ind);
					    fitness = resizeArray(fitness, worst_ind);
					    LinearPopSize--;
					    //arc_size--;
					}
					arc_size = (int)Math.round(LinearPopSize * arc_rate);
					//Psxhll: update bestFitness value and its index, as we change the size after finishing each for loop
					bestFitness = fitness[0];
					bestFitnessIndex =0;
					for(int i=1; i< population.length;i++){
						  if(fitness[i]>bestFitness){
							  bestFitness = fitness[i];
							  bestFitnessIndex=i;
						  }
						  
					 }
					
			    }
	   } // End main LOOP
	     
	//   System.out.println("training Size" +trainingDataSet.size());
	   System.out.println("Best Fitness "+bestFitness);
		   nominalPopulation = new PrototypeSet();
           nominalPopulation.formatear(population[bestFitnessIndex]);
           
  //         System.out.println("Best Fitness2 "+  accuracy(nominalPopulation,trainingDataSet));
			 System.err.println("\n% de acierto en training Nominal " + KNN.classficationAccuracy(nominalPopulation,trainingDataSet,1)*100./trainingDataSet.size() );
				  
			//  nominalPopulation.print();

  
		return nominalPopulation;
  }
  
  
  public double[] resizeArray(double[] ori, int del_index) {
		double[] resized_array = new double[ori.length - 1];
		int index_count = 0;

		for (int i = 0; i < ori.length; i++) {
		    if (i != del_index) {
		    	resized_array[index_count] = ori[i];
		    	index_count++;
		    }
		}

		return resized_array;
  }

 public PrototypeSet[] resizeMatrix(PrototypeSet[] ori, int del_index) {
	   PrototypeSet[] newset = new PrototypeSet[ori.length-1];
	   int index_count = 0;
	   for (int i = 0; i < ori.length; i++) {
		    if (i != del_index) {
		    	newset[index_count] = ori[i].clone();
		    	index_count++;
		    }
		}
	   return newset;
 }
  
  /**
   * General main for all the prototoype generators
   * Arguments:
   * 0: Filename with the training data set to be condensed.
   * 1: Filename which contains the test data set.
   * 3: Seed of the random number generator.            Always.
   * **************************
   * 4: .Number of neighbors
   * 5:  Swarm Size
   * 6:  Particle Size
   * 7:  Max Iter
   * 8:  C1
   * 9: c2
   * 10: vmax
   * 11: wstart
   * 12: wend
   * @param args Arguments of the main function.
 * @throws Exception 
   */
  public static void main(String[] args) throws Exception
  {
      Parameters.setUse("JADE", "<seed> <Number of neighbors>\n<Swarm size>\n<Particle Size>\n<MaxIter>\n<DistanceFunction>");        
      Parameters.assertBasicArgs(args);
      
      PrototypeSet training = PrototypeGenerationAlgorithm.readPrototypeSet(args[0]);
      PrototypeSet test = PrototypeGenerationAlgorithm.readPrototypeSet(args[1]);
      
      
      long seed = Parameters.assertExtendedArgAsInt(args,2,"seed",0,Long.MAX_VALUE);
      LSHADEGenerator.setSeed(seed);
      
      int k = Parameters.assertExtendedArgAsInt(args,3,"number of neighbors", 1, Integer.MAX_VALUE);
      int swarm = Parameters.assertExtendedArgAsInt(args,4,"swarm size", 1, Integer.MAX_VALUE);
      int particle = Parameters.assertExtendedArgAsInt(args,5,"particle size", 1, Integer.MAX_VALUE);
      int iter = Parameters.assertExtendedArgAsInt(args,6,"max iter", 1, Integer.MAX_VALUE);
      double c1 = Parameters.assertExtendedArgAsInt(args,7,"c1", 1, Double.MAX_VALUE);
      double c2 =Parameters.assertExtendedArgAsInt(args,8,"c2", 1, Double.MAX_VALUE);
      double vmax =Parameters.assertExtendedArgAsInt(args,9,"vmax", 1, Double.MAX_VALUE);
      double wstart = Parameters.assertExtendedArgAsInt(args,10,"wstart", 1, Double.MAX_VALUE);
      double wend =Parameters.assertExtendedArgAsInt(args,11,"wend", 1, Double.MAX_VALUE);
      
      //String[] parametersOfInitialReduction = Arrays.copyOfRange(args, 4, args.length);
     //System.out.print(" swarm ="+swarm+"\n");
      
      
      LSHADEGenerator generator = new LSHADEGenerator(training, k,swarm,particle,iter, 0.5,0.5,1);
      
  	  
      PrototypeSet resultingSet = generator.execute();
      
  	//resultingSet.save(args[1]);
      //int accuracyKNN = KNN.classficationAccuracy(resultingSet, test, k);
      int accuracy1NN = KNN.classficationAccuracy(resultingSet, test);
      generator.showResultsOfAccuracy(Parameters.getFileName(), accuracy1NN, test);
  }

}
