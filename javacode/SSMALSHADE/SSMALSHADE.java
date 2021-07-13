//
//  SSMA.java
//
//  Salvador Garcï¿½a Lï¿½pez
//
//  Created by Salvador Garcï¿½a Lï¿½pez 3-10-2005.
//  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
//

package keel.Algorithms.Instance_Generation.SSMALSHADE;

import keel.Algorithms.Preprocess.Basic.*;

import keel.Algorithms.Instance_Generation.Basic.Prototype;
import keel.Algorithms.Instance_Generation.Basic.PrototypeGenerationAlgorithm;
import keel.Algorithms.Instance_Generation.Basic.PrototypeSet;
import keel.Algorithms.Instance_Generation.Basic.PublicVar;
import keel.Algorithms.Instance_Generation.utilities.*;
import keel.Algorithms.Instance_Generation.Basic.PrototypeGenerator;



import keel.Dataset.Attribute;
import keel.Dataset.Attributes;
import keel.Dataset.InstanceAttributes;
import keel.Dataset.InstanceSet;

import org.core.*;

import java.util.StringTokenizer;
import java.util.Arrays;

public class SSMALSHADE extends Metodo {

  /*Own parameters of the algorithm*/

  private long semilla;
  private int tamPoblacion;
  private double nEvalSSMA;
  private double pCross;
  private double pMut;
  private int kNeigh;
  public String Script; // para releer parï¿½metros..
  private PrototypeSet trainingDataSet;
  private PrototypeSet testDataSet;
  private PrototypeGenerator generador;
  protected int numberOfbetters; // numero de mejores atener en cuenta
  //Parï¿½metros DE
  private int k;
  
  private int PopulationSize; 
  private int ParticleSize;
  private int maxEval_IG; 
  private double ScalingFactor;
  private double CrossOverRate[];
  private int Strategy;
  private String CrossoverType; // Binomial, Exponential, Arithmetic
  
  
  private double tau[] = new double[4];
  private double Fl, Fu;
  
  private int iterSFGSS;
  private int iterSFHC;
  
  protected int numberOfClass;


  protected int numberOfPrototypes;  // Particle size is the percentage
  protected int numberOfStrategies; // number of strategies in the pool
  private int perc;
  

  public SSMALSHADE (String ficheroScript) {
	super (ficheroScript);
    
  }


  public SSMALSHADE(String ficheroScript, InstanceSet train) {
		super (ficheroScript, train);
  }

 public void establishTrain(PrototypeSet trainPG){
	 trainingDataSet = trainPG.clone();
 }

/**
   * Reads the prototype set from a data file.
   * @param nameOfFile Name of data file to be read.
   * @return PrototypeSet built with the data of the file.
   */
  public static PrototypeSet readPrototypeSet(String nameOfFile)
  {
      Attributes.clearAll();
      InstanceSet training = new InstanceSet();        
      try
      {
      	//System.out.print("PROBANDO:\n"+nameOfFile);
          training.readSet(nameOfFile, true); 
          training.setAttributesAsNonStatic();
          InstanceAttributes att = training.getAttributeDefinitions();
          Prototype.setAttributesTypes(att);            
      }
      catch(Exception e)
      {
          System.err.println("readPrototypeSet has failed!");
          e.printStackTrace();
      }
      return new PrototypeSet(training);
  }
  
  
  public static PrototypeSet readPrototypeSet2(InstanceSet training)
  {
      Attributes.clearAll();
 
      try
      {
          training.setAttributesAsNonStatic();
          InstanceAttributes att = training.getAttributeDefinitions();
          Prototype.setAttributesTypes(att);            
      }
      catch(Exception e)
      {
          System.err.println("readPrototypeSet has failed!");
          e.printStackTrace();
      }
      return new PrototypeSet(training);
  }
  
  
  
  public void inic_vector_sin(int vector[], int without){

  	for(int i=0; i<vector.length; i++) 
  		if(i!=without)
  			vector[i] = i; 
  }
  
  public void desordenar_vector_sin(int vector[]){
  	int tmp, pos;
  	for(int i=0; i<vector.length-1; i++){
  		pos = Randomize.Randint(0, vector.length-1);
  		tmp = vector[i];
  		vector[i] = vector[pos];
  		vector[pos] = tmp;
  	}
  }
  

  
  public PrototypeSet mutant(PrototypeSet population[], int actual, int mejor, double SFi){
  	  
  	  
  	  PrototypeSet mutant = new PrototypeSet(population.length);
  	  PrototypeSet r1,r2,r3,r4,r5, resta, producto, resta2, producto2, result, producto3, resta3;
  	  
  	//We need three differents solutions of actual
  		   
  	  int lista[] = new int[population.length];
        inic_vector_sin(lista,actual);
        desordenar_vector_sin(lista);
  		      
  	  // System.out.println("Lista = "+lista[0]+","+ lista[1]+","+lista[2]);
  	  
  	   r1 = population[lista[0]];
  	   r2 = population[lista[1]];
  	   r3 = population[lista[2]];
  	   r4 = population[lista[3]];
  	   r5 = population[lista[4]];
  		   
  			switch(this.Strategy){
  		   	   case 1: // ViG = Xr1,G + F(Xr2,G - Xr3,G) De rand 1
  		   		 resta = r2.restar(r3); //Psxhll: r2 - r3
  		   		 producto = resta.mulEscalar(SFi);
  		   		 mutant = producto.sumar(r1);
  		   	    break;
  			   
  		   	   case 2: // Vig = Xbest,G + F(Xr2,G - Xr3,G)  De best 1
  			   		 resta = r2.restar(r3);
  			   		 producto = resta.mulEscalar(SFi);
  			   		 mutant = population[mejor].sumar(producto);
  			   break;
  			   
  		   	   case 3: // Vig = ... De rand to best 1
  		   		   resta = r1.restar(r2); 
  		   		   resta2 = population[mejor].restar(population[actual]);
  		   		 			   		 
  			   	   producto = resta.mulEscalar(SFi);
  			   	   producto2 = resta2.mulEscalar(SFi);
  			   		
  			   	   result = population[actual].sumar(producto);
  			   	   mutant = result.sumar(producto2);
  			   		 			   		 
  			   break;
  			   
  		   	   case 4: // DE best 2
  		   		   resta = r1.restar(r2); 
  		   		   resta2 = r3.restar(r4);
  		   		 			   		 
  			   	   producto = resta.mulEscalar(SFi);
  			   	   producto2 = resta2.mulEscalar(SFi);
  			   		
  			   	   result = population[mejor].sumar(producto);
  			   	   mutant = result.sumar(producto2);
  			   break;
  			  
  		   	   case 5: //DE rand 2
  		   		   resta = r2.restar(r3); 
  		   		   resta2 = r4.restar(r5);
  		   		 			   		 
  			   	   producto = resta.mulEscalar(SFi);
  			   	   producto2 = resta2.mulEscalar(SFi);
  			   		
  			   	   result = r1.sumar(producto);
  			   	   mutant = result.sumar(producto2);
  			   	   
    		       break;
    		       
  		   	   case 6: //DE rand to best 2
  		   		   resta = r1.restar(r2); 
  		   		   resta2 = r3.restar(r4);
  		   		   resta3 = population[mejor].restar(population[actual]);
  		   		   
  			   	   producto = resta.mulEscalar(SFi);
  			   	   producto2 = resta2.mulEscalar(SFi);
  			   	   producto3 = resta3.mulEscalar(SFi);
  			   	   
  			   	   result = population[actual].sumar(producto);
  			   	   result = result.sumar(producto2);
  			   	   mutant = result.sumar(producto3);
    		       break;
    		       
  		   }   
  	   

  	  // System.out.println("********Mutante**********");
  	 // mutant.print();
  	   
       mutant.applyThresholds();
  	
  	  return mutant;
    }


  /**
   * Local Search Fitness Function
   * @param Fi
   * @param xt
   * @param xr
   * @param xs
   * @param actual
   */
  public double lsff(double Fi, double CRi, PrototypeSet population[], int actual, int mejor){
	  PrototypeSet resta, producto, mutant;
	  PrototypeSet crossover;
	  double FitnessFi = 0;
	  
	  
	  //Mutation:
	  mutant = new PrototypeSet(population[actual].size());
   	  mutant = mutant(population, actual, mejor, Fi);
   	
   	  
   	  //Crossover
   	  crossover =new PrototypeSet(population[actual]);
   	  
	   for(int j=0; j< population[actual].size(); j++){ // For each part of the solution
		   
		   double randNumber = RandomGenerator.Randdouble(0, 1);
			   
		   if(randNumber< CRi){
			   crossover.set(j, mutant.get(j)); // Overwrite.
		   }
	   }
	   
	   
	   // Compute fitness
	   PrototypeSet nominalPopulation = new PrototypeSet();
           nominalPopulation.formatear(crossover);
           FitnessFi =  classficationAccuracy1NN(nominalPopulation,trainingDataSet);
	   
   	   return FitnessFi;
  }
  
  
  
  /**
   * SFGSS local Search.
   * @param population
   * @return
   */
  public PrototypeSet SFGSS(PrototypeSet population[], int actual, int mejor, double CRi){
	  double a=0.1, b=1;

	  double fi1=0, fi2=0, fitnessFi1=0, fitnessFi2=0;
	  double phi = (1+ Math.sqrt(5))/5;
	  double scaling;
	  PrototypeSet crossover, resta, producto, mutant;
	  
	  for (int i=0; i<this.iterSFGSS; i++){ // Computation budjet
	  
		  fi1 = b - (b-a)/phi;
		  fi2 = a + (b-a)/phi;
		  
		  fitnessFi1 = lsff(fi1, CRi, population,actual,mejor);
		  fitnessFi2 = lsff(fi2, CRi,population,actual,mejor);
		  
		  if(fitnessFi1> fitnessFi2){
			  b = fi2;
		  }else{
			  a = fi1;  
		  }
	  
	  } // End While
	  
	  
	  if(fitnessFi1> fitnessFi2){
		  scaling = fi1;
	  }else{
		  scaling = fi2;
	  }
	  
	  
	  //Mutation:
	  mutant = new PrototypeSet(population[actual].size());
	  mutant = mutant(population, actual, mejor, scaling);
   	  
   	  //Crossover
   	  crossover =new PrototypeSet(population[actual]);
   	  
	   for(int j=0; j< population[actual].size(); j++){ // For each part of the solution
		   
		   double randNumber = RandomGenerator.Randdouble(0, 1);
			   
		   if(randNumber< CRi){
			   crossover.set(j, mutant.get(j)); // Overwrite.
		   }
	   }
	   
	   
	  
	return crossover;
  }
  
  /**
   * SFHC local search
   * @param xt
   * @param xr
   * @param xs
   * @param actual
   * @param SFi
   * @return
   */
  
  public  PrototypeSet SFHC(PrototypeSet population[], int actual, int mejor, double SFi, double CRi){
	  double fitnessFi1, fitnessFi2, fitnessFi3, bestFi;
	  PrototypeSet crossover, resta, producto, mutant;
	  double h= 0.5;
	  
	  
	  for (int i=0; i<this.iterSFHC; i++){ // Computation budjet
		  		  
		  fitnessFi1 = lsff(SFi-h, CRi, population,actual,mejor);
		  fitnessFi2 = lsff(SFi, CRi,  population,actual,mejor);
		  fitnessFi3 = lsff(SFi+h, CRi,  population,actual,mejor);
		  
		  if(fitnessFi1 >= fitnessFi2 && fitnessFi1 >= fitnessFi3){
			  bestFi = SFi-h;
		  }else if(fitnessFi2 >= fitnessFi1 && fitnessFi2 >= fitnessFi3){
			  bestFi = SFi;
			  h = h/2; // H is halved.
		  }else{
			  bestFi = SFi;
		  }
		  
		  SFi = bestFi;
	  }
	  
	  
	  //Mutation:
	  mutant = new PrototypeSet(population[actual].size());
	  mutant = mutant(population, actual, mejor, SFi);
	 
   	  //Crossover
   	  crossover = new PrototypeSet(population[actual]);
   	  
	   for(int j=0; j< population[actual].size(); j++){ // For each part of the solution
		   
		   double randNumber = RandomGenerator.Randdouble(0, 1);
			   
		   if(randNumber< CRi){
			   crossover.set(j, mutant.get(j)); // Overwrite.
		   }
	   }
	   
	   
	  
	return crossover;
  
  }
  
  
  
  /**
   * Implements the 1NN algorithm
   * @param current Prototype which the algorithm will find its nearest-neighbor.
   * @param dataSet Prototype set in which the algorithm will search.
   * @return Nearest prototype to current in the prototype set dataset.
   */
  public static Prototype _1nn(Prototype current, PrototypeSet dataSet)
  {
      Prototype nearestNeighbor = dataSet.get(0);
      int indexNN = 0;
      //double minDist = Distance.dSquared(current, nearestNeighbor);
      //double minDist = Distance.euclideanDistance(current, nearestNeighbor);
      double minDist =Double.POSITIVE_INFINITY;
      double currDist;
      int _size = dataSet.size();
    //  System.out.println("****************");
     // current.print();
      for (int i=0; i<_size; i++)
      {
          Prototype pi = dataSet.get(i);
          //if(!current.equals(pi))
          //{
             // double currDist = Distance.dSquared(current, pi);
           currDist = Distance.euclideanDistance(pi,current);
          // System.out.println(currDist);
          
           if(currDist >0){
              if (currDist < minDist)
              {
                  minDist = currDist;
                 // nearestNeighbor = pi;
                  indexNN =i;
              }
          }
          //}
      }
      
     // System.out.println("Min dist =" + minDist + " Vecino Cercano = "+ indexNN);
      
      return dataSet.get(indexNN);
  }
  


//  Sarah: TO BE MODIFIED!!!!!

  public double classficationAccuracy1NN(PrototypeSet training, PrototypeSet test)
  {
	int wellClassificated = 0;
      for(Prototype p : test)
      {
          Prototype nearestNeighbor = _1nn(p, training);          
          
          if(p.getOutput(0) == nearestNeighbor.getOutput(0))
              ++wellClassificated;
      }
  
      
      return 100.0* (wellClassificated / (double)test.size());
  }
  
  
  
  
  /**
   * Generate a reduced prototype set by the SADEGenerator method.
   * @return Reduced set by SADEGenerator's method.
   */
  

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
	  
	  
	  public PrototypeSet reduceSet(PrototypeSet initial)
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
		  population[0]= new PrototypeSet(initial.clone()) ;
		  		  
		  generador = new PrototypeGenerator(trainingDataSet);
		  
		  // If SSMA fails (It may happen, I know!):
		  if(population[0].size() <2){
			  this.numberOfPrototypes = (int)Math.round(trainingDataSet.size() * this.perc/100);
			  population[0]=generador.selecRandomSet(numberOfPrototypes,true).clone() ;
			  // red .95
			  
			  // Aseguro que al menos hay un representante de cada clase.
			  PrototypeSet clases[] = new PrototypeSet [this.numberOfClass];
			  for(int i=0; i< this.numberOfClass; i++){
				  clases[i] = new PrototypeSet(trainingDataSet.getFromClass(i));
				  
				 // System.out.println("Clase "+i+", size= "+ clases[i].size());
			  }
			
			  for(int i=0; i< population[0].size(); i++){
				  for(int j=0; j< this.numberOfClass; j++){
					  if(population[0].getFromClass(j).size() ==0 && clases[j].size()!=0){
						  
						  population[0].add(clases[j].getRandom());
					  }
				  }
			  }
		  }
		  fitness[0] = classficationAccuracy1NN(population[0],trainingDataSet);
		  
		 // population[0].print();
		  for(int i=1; i< PopulationSize; i++){
			  population[i] = new PrototypeSet();
			  for(int j=0; j< population[0].size(); j++){
				  population[i].add(trainingDataSet.getFromClass(population[0].get(j).getOutput(0)).getRandom());
			  }
			  fitness[i] = classficationAccuracy1NN(population[i],trainingDataSet);   // DE fitness, no hace falta formatear porque son aleatorios!
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
		   while (ev < maxEval_IG)// Main loop
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
				       fitness[i] = classficationAccuracy1NN(nominalPopulation,trainingDataSet);
				       ev++;	  
				       nominalPopulation = new PrototypeSet();
				       nominalPopulation.formatear(crossover[i]);
				       
					  double trialVector = classficationAccuracy1NN(nominalPopulation,trainingDataSet);
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
				  int next_pop_size = (int)Math.round((((min_pop_size - PopulationSize) / (double)maxEval_IG) * ev) + PopulationSize);
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
				 System.err.println("\n% de acierto en training Nominal " + classficationAccuracy1NN(nominalPopulation,trainingDataSet));
					  
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
	  
  public PrototypeSet reduceSet_OriginalSetting(PrototypeSet initial)
  {
	  	  
	  System.out.print("\nThe algorithm  SSMA-SFLSDE is starting...\n Computing...\n");
	  
	  //Algorithm
	  // First, we create the population, with PopulationSize.
	  // like a prototypeSet's vector.  

	  PrototypeSet population [] = new PrototypeSet [PopulationSize];
	  PrototypeSet mutation[] = new PrototypeSet[PopulationSize];
	  PrototypeSet crossover[] = new PrototypeSet[PopulationSize];
	  
	  
	  double ScalingFactor[] = new double[this.PopulationSize];
	  double CrossOverRate[] = new double[this.PopulationSize]; // Inside of the Optimization process.
	  double fitness[] = new double[PopulationSize];

	  double fitness_bestPopulation[] = new double[PopulationSize];
	  PrototypeSet bestParticle = new PrototypeSet();
	
  
	  //Each particle must have   Particle Size %

	  // First Stage, Initialization.
	  
	  PrototypeSet nominalPopulation;
	  
	  population[0]= new PrototypeSet(initial.clone()) ;
	  
	  generador = new PrototypeGenerator(trainingDataSet);
	  
	  
	  // If SSMA fails (It may happen, I know!):
	  
	 // population[0].print();
	  
	  if(population[0].size() <2){
		  this.numberOfPrototypes = (int)Math.round(trainingDataSet.size()*0.02);
		  
		  population[0]=generador.selecRandomSet(numberOfPrototypes,true).clone() ;
		  // red .95
		  
		  // Aseguro que al menos hay un representante de cada clase.
		  PrototypeSet clases[] = new PrototypeSet [this.numberOfClass];
		  for(int i=0; i< this.numberOfClass; i++){
			  clases[i] = new PrototypeSet(trainingDataSet.getFromClass(i));
			  
			 // System.out.println("Clase "+i+", size= "+ clases[i].size());
		  }
		
		  for(int i=0; i< population[0].size(); i++){
			  for(int j=0; j< this.numberOfClass; j++){
				  if(population[0].getFromClass(j).size() ==0 && clases[j].size()!=0){
					  
					  population[0].add(clases[j].getRandom());
				  }
			  }
		  }
	  }
	  
	  
	
	  nominalPopulation = new PrototypeSet();
          nominalPopulation.formatear(population[0]);
       
	  fitness[0] = classficationAccuracy1NN(nominalPopulation,trainingDataSet);
	  
	  System.out.println("Best initial fitness = "+ fitness[0]);

	  this.numberOfClass = trainingDataSet.getPosibleValuesOfOutput().size();
	  
	  
      for(int i=1; i< PopulationSize; i++){
		  population[i] = new PrototypeSet();
		 
		  for(int j=0; j< population[0].size(); j++){
			  Prototype aux = new Prototype(trainingDataSet.getFromClass(population[0].get(j).getOutput(0)).getRandom());
			  population[i].add(aux);
		  }
		  
		  nominalPopulation = new PrototypeSet();
	          nominalPopulation.formatear(population[i]);
	      
		  fitness[i] = classficationAccuracy1NN(population[i],trainingDataSet);   // PSOfitness
		  fitness_bestPopulation[i] = fitness[i]; // Initially the same fitness.
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
	   
	   
	   // Initially the Scaling Factor and crossover for each Individual are randomly generated between 0 and 1.
	   
	   for(int i=0; i< this.PopulationSize; i++){
		   ScalingFactor[i] =  RandomGenerator.Randdouble(0, 1);
		   CrossOverRate[i] =  RandomGenerator.Randdouble(0, 1);
	   }
	   
	   double randj[] = new double[5];
	   
	   
	   for(int iter=0; iter< maxEval_IG; iter++){ // Main loop
	   
		   for(int i=0; i<PopulationSize; i++)
		   {
			   // Generate randj for j=1 to 5.
			   for(int j=0; j<5; j++){
				   randj[j] = RandomGenerator.Randdouble(0, 1);
			   }
			   
	   
			   if(i==bestFitnessIndex && randj[4] < tau[2]){
				   System.out.println("SFGSS applied");
				   //SFGSS
				   crossover[i] = SFGSS(population, i, bestFitnessIndex, CrossOverRate[i]);
				   
			   }else if(i==bestFitnessIndex &&  tau[2] <= randj[4] && randj[4] < tau[3]){
				   //SFHC
				   System.out.println("SFHC applied");
				   crossover[i] = SFHC(population, i, bestFitnessIndex, ScalingFactor[i], CrossOverRate[i]);
				   
			   }else {
				   //System.out.println("Random applied");
				   // Fi update
				   
				   if(randj[1] < tau[0]){
					   ScalingFactor[i] = this.Fl + this.Fu*randj[0];
				   }
				   
				   // CRi update
				   
				   if(randj[3] < tau[1]){
					   CrossOverRate[i] = randj[2];
				   }
				   				   
				   // Mutation Operation.
				   
				   mutation[i] = new PrototypeSet(population[i].size());
			   
				  //Mutation:
					
				   mutation[i]  = mutant(population, i, bestFitnessIndex, ScalingFactor[i]);
				   
				    // Crossver Operation.

				   crossover[i] = new PrototypeSet(population[i]);
				   
				   for(int j=0; j< population[i].size(); j++){ // For each part of the solution
					   
					   double randNumber = RandomGenerator.Randdouble(0, 1);
						   
					   if(randNumber<CrossOverRate[i]){
						   crossover[i].set(j, mutation[i].get(j)); // Overwrite.
					   }
				   }
								   
			   }
				// Fourth: Selection Operation.
		   
			   // nominalPopulation = new PrototypeSet();
		           // nominalPopulation.formatear(population[i]);
		           // fitness[i] = classficationAccuracy1NN(nominalPopulation,trainingDataSet);
		       
	            nominalPopulation = new PrototypeSet();
	            nominalPopulation.formatear(crossover[i]);
		       
			    double trialVector = classficationAccuracy1NN(nominalPopulation,trainingDataSet);
		  
			  if(trialVector > fitness[i]){
				  population[i] = new PrototypeSet(crossover[i]);
				  fitness[i] = trialVector;
			  }
			  
			  if(fitness[i]>bestFitness){
				  bestFitness = fitness[i];
				  bestFitnessIndex=i;
			  }
			  
			  
		   }

		   
	   }//End main while-loop

	   nominalPopulation = new PrototypeSet();
           nominalPopulation.formatear(population[bestFitnessIndex]);
	   System.err.println("\n% de acierto en training Nominal " + classficationAccuracy1NN(nominalPopulation,trainingDataSet) );
			  
	   //  nominalPopulation.print();
	   return nominalPopulation;
  }
  
  
  /* MEzcla de algoritmos */
  public void ejecutar () {
	double count =0, countLocal=0;
    int i, j, l;
    double conjS[][];
    double conjR[][];
    int conjN[][];
    boolean conjM[][];
    int clasesS[];
    int nSel = 0;
    Cromosoma poblacion[];
    double ev = 0;
    double dMatrix[][];
    int sel1, sel2, comp1, comp2;
    Cromosoma hijos[];
    double umbralOpt;
    boolean veryLarge;
    double GAeffort=0, LSeffort=0, temporal;
    double fAcierto=0, fReduccion=0;
    int contAcierto=0, contReduccion=0;
    int nClases;

    long tiempo = System.currentTimeMillis();

    /*Getting the number of different classes*/
    nClases = 0;
    for (i=0; i<clasesTrain.length; i++)
      if (clasesTrain[i] > nClases)
        nClases = clasesTrain[i];
    nClases++;

    if (datosTrain.length > 9000) {
      veryLarge = true;
    } else {
      veryLarge = false;
    }

    if (veryLarge == false) {
      /*Construct a distance matrix of the instances*/
      dMatrix = new double[datosTrain.length][datosTrain.length];
      for (i = 0; i < dMatrix.length; i++) {
        for (j = i + 1; j < dMatrix[i].length; j++) {
          dMatrix[i][j] = KNN.distancia(datosTrain[i], realTrain[i], nominalTrain[i], nulosTrain[i], datosTrain[j], realTrain[j], nominalTrain[j], nulosTrain[j], distanceEu);
        }
      }
      for (i = 0; i < dMatrix.length; i++) {
        dMatrix[i][i] = Double.POSITIVE_INFINITY;
      }
      for (i = 0; i < dMatrix.length; i++) {
        for (j = i - 1; j >= 0; j--) {
          dMatrix[i][j] = dMatrix[j][i];
        }
      }
    } else {
      dMatrix = null;
    }

    /*Random inicialization of the population*/
    Randomize.setSeed (semilla);
    poblacion = new Cromosoma[tamPoblacion];
    for (i=0; i<tamPoblacion; i++)
      poblacion[i] = new Cromosoma (kNeigh, datosTrain.length, dMatrix, datosTrain, realTrain, nominalTrain, nulosTrain, distanceEu);
	int activegenBefore = poblacion[0].genesActivos();
    /*Initial evaluation of the population*/
    for (i=0; i<tamPoblacion; i++) {
      poblacion[i].evaluacionCompleta(nClases, kNeigh, clasesTrain);
    }

    umbralOpt = 0;

    /*Until stop condition*/
    //Psxhll note: As SSMA is used in a combination with other instance generation techniques. The number of evaluation is modified.
    if (PublicVar.OriginalSetting)
    	nEvalSSMA = PublicVar.OriginalEvalSSMA;
    else
    	nEvalSSMA = PublicVar.EvalWeight * datosTrain.length;
    	
    while (ev < nEvalSSMA) {

      Arrays.sort(poblacion);

      if (fAcierto >= (double)poblacion[0].getFitnessAc()*100.0/(double)datosTrain.length) {
        contAcierto++;
      } else {
        contAcierto=0;
      }
      fAcierto = (double)poblacion[0].getFitnessAc()*100.0/(double)datosTrain.length;

      if (fReduccion >= (1.0-((double)poblacion[0].genesActivos()/(double)datosTrain.length))*100.0) {
        contReduccion++;
      } else {
        contReduccion=0;
      }
      fReduccion = (1.0-((double)poblacion[0].genesActivos()/(double)datosTrain.length))*100.0;

      if (contReduccion >= 10 || contAcierto >= 10){
        if (Randomize.Randint(0,1)==0) {
          if (contAcierto >= 10) {
            contAcierto = 0;
            umbralOpt++;
          } else {
            contReduccion = 0;
            umbralOpt--;
          }
        } else {
          if (contReduccion >= 10) {
            contReduccion = 0;
            umbralOpt--;
          } else {
            contAcierto = 0;
            umbralOpt++;
          }
        }
      }

      /*Binary tournament selection*/
      comp1 = Randomize.Randint(0,tamPoblacion-1);
      do {
        comp2 = Randomize.Randint(0,tamPoblacion-1);
      } while (comp2 == comp1);

      if (poblacion[comp1].getFitness() > poblacion[comp2].getFitness())
        sel1 = comp1;
      else sel1 = comp2;
      comp1 = Randomize.Randint(0,tamPoblacion-1);
      do {
        comp2 = Randomize.Randint(0,tamPoblacion-1);
      } while (comp2 == comp1);
      if (poblacion[comp1].getFitness() > poblacion[comp2].getFitness())
        sel2 = comp1;
      else
        sel2 = comp2;


      hijos = new Cromosoma[2];
      hijos[0] = new Cromosoma (kNeigh, poblacion[sel1], poblacion[sel2], pCross,datosTrain.length);
      hijos[1] = new Cromosoma (kNeigh, poblacion[sel2], poblacion[sel1], pCross,datosTrain.length);
      hijos[0].mutation (kNeigh, pMut, dMatrix, datosTrain, realTrain, nominalTrain, nulosTrain, distanceEu);
      hijos[1].mutation (kNeigh, pMut, dMatrix, datosTrain, realTrain, nominalTrain, nulosTrain, distanceEu);

      /*Evaluation of offsprings*/
	      long time_begin = System.currentTimeMillis();
                    
	      hijos[0].evaluacionCompleta(nClases, kNeigh, clasesTrain);
	      hijos[1].evaluacionCompleta(nClases, kNeigh, clasesTrain);
	      double time_end = (double) (System.currentTimeMillis() - time_begin) / 1000.0;
	      count +=2;
	      ev+=2;
	      GAeffort += 2;
	      temporal = ev;
	      if (hijos[0].getFitness() > poblacion[tamPoblacion-1].getFitness() || Randomize.Rand() < 0.0625) {
	    	  double localev =  hijos[0].optimizacionLocal(nClases, kNeigh, clasesTrain,dMatrix,umbralOpt, datosTrain, realTrain, nominalTrain, nulosTrain, distanceEu);
	    	  ev += localev;
	    	  countLocal += localev;
	      }
	      
	      if (hijos[1].getFitness() > poblacion[tamPoblacion-1].getFitness() || Randomize.Rand() < 0.0625) {
	    	  double localev = hijos[1].optimizacionLocal(nClases, kNeigh, clasesTrain,dMatrix,umbralOpt, datosTrain, realTrain, nominalTrain, nulosTrain, distanceEu);
	          ev += localev;
	    	  countLocal += localev;
	      }
	      double time_end2 = (double) (System.currentTimeMillis() - time_begin) / 1000.0;
	      fitCollect.timeComplete_Local.add(new Pair<Double, Double> (time_end, time_end2));
	      LSeffort += (ev - temporal);

      /*Replace the two worst*/
      if (hijos[0].getFitness() > poblacion[tamPoblacion-1].getFitness()) {
        poblacion[tamPoblacion-1] = new Cromosoma (kNeigh, datosTrain.length, hijos[0]);
      }
      if (hijos[1].getFitness() > poblacion[tamPoblacion-2].getFitness()) {
        poblacion[tamPoblacion-2] = new Cromosoma (kNeigh, datosTrain.length, hijos[1]);
      }
      /*System.out.println(ev + " - (" + umbralOpt + ")" +
                         (double)poblacion[0].getFitnessAc()*100.0/(double)datosTrain.length + " / " +
                         (double)poblacion[tamPoblacion-1].getFitnessAc()*100.0/(double)datosTrain.length + " - " +
                         (1.0-((double)poblacion[0].genesActivos()/(double)datosTrain.length))*100.0 + " / " +
                         (1.0-((double)poblacion[tamPoblacion-1].genesActivos()/(double)datosTrain.length))*100.0);*/
    }//End SSMA

    Arrays.sort(poblacion);
    nSel = poblacion[0].genesActivos();

    /*Construction of S set from the best cromosome*/
    conjS = new double[nSel][datosTrain[0].length];
    conjR = new double[nSel][datosTrain[0].length];
    conjN = new int[nSel][datosTrain[0].length];
    conjM = new boolean[nSel][datosTrain[0].length];
    clasesS = new int[nSel];
    for (i=0, l=0; i<datosTrain.length; i++) {
      if (poblacion[0].getGen(i)) { //the instance must be copied to the solution
        for (j=0; j<datosTrain[i].length; j++) {
          conjS[l][j] = datosTrain[i][j];
          conjR[l][j] = realTrain[i][j];
          conjN[l][j] = nominalTrain[i][j];
          conjM[l][j] = nulosTrain[i][j];
        }
        clasesS[l] = clasesTrain[i];
        l++;
      }
    }

    System.out.println("SSMA has done " + (int)(ev) + " evaluations. Max evaluation: " + nEvalSSMA +"\n\n");
	double ssmaruntime = (double)(System.currentTimeMillis()-tiempo)/1000.0;
	System.out.println("--- Runtime of SSMA  for "+ relation + " " + ssmaruntime + "s");

    OutputIS.escribeSalida(ficheroSalida[0], conjR, conjN, conjM, clasesS, entradas, salida, nEntradas, relation);
    OutputIS.escribeSalida(ficheroSalida[1], test, entradas, salida, nEntradas, relation);
    
    

    Parameters.assertBasicArgs(ficheroSalida);
    
    if(!this.Script.equals("NOFILE")){
    	PrototypeGenerationAlgorithm.readParametersFile(this.Script);
    	PrototypeGenerationAlgorithm.printParameters();
        trainingDataSet = readPrototypeSet(this.ficheroTraining); 
        testDataSet = readPrototypeSet(this.ficheroTest);
    }

    PrototypeSet training = readPrototypeSet(ficheroSalida[0]);
    //--------------------------------------------------------
	// If SSMA fails (It may happen, I know!):
	 	  
	  if(training.size() <2){
		  this.numberOfPrototypes = (int)Math.round(trainingDataSet.size() * this.perc/100);
		  generador = new PrototypeGenerator(trainingDataSet);
		  training = generador.selecRandomSet(numberOfPrototypes, true).clone() ;
		  // red .95
		  
		  // Aseguro que al menos hay un representante de cada clase.
		  this.numberOfClass = trainingDataSet.getPosibleValuesOfOutput().size();
		  PrototypeSet clases[] = new PrototypeSet [this.numberOfClass];
		  for(i=0; i< this.numberOfClass; i++){
			  clases[i] = new PrototypeSet(trainingDataSet.getFromClass(i));
			  
			 // System.out.println(""Clase ""+i+"", size= ""+ clases[i].size());
		  }
		
		  for(i=0; i< training.size(); i++){
			  for(j=0; j< this.numberOfClass; j++){
				  if(training.getFromClass(j).size() ==0 && clases[j].size()!=0){
					  
					  training.add(clases[j].getRandom());
				  }
			  }
		  }
	  }
	  //--------------------------------------------------------"

   // trainingDataSet.print();
   //this.numberOfPrototypes = (int)Math.floor((trainingDataSet.size())*ParticleSize/100.0);
   // training.print();
    PrototypeSet SADE;
    if (PublicVar.OriginalSetting) {
    	this.maxEval_IG = PublicVar.OriginalSFLSDEMaxIter;
    	SADE = reduceSet_OriginalSetting(training);
    }
    else {
    	int dimensions = training.get(0).numberOfInputs();
 	    this.maxEval_IG = (int) (trainingDataSet.size() * (dimensions*PublicVar.EvalWeight - PublicVar.SSMAWeight));
    	SADE = reduceSet (training);
    }
	//Psxhll: Print out the fitness collected in child class. The child class needs to override the method in the parent class
	    fitCollect.writeToFile(ficheroSalida[0]);
	    fitCollect.writeToFile_ReductionRate(ficheroSalida[0], training.size(), trainingDataSet.size());
	    double IGruntime = (double)(System.currentTimeMillis()-tiempo)/1000.0 - ssmaruntime;   
	    fitCollect.writeToFile_Time_SSMA_IG(ficheroSalida[0], ssmaruntime, IGruntime, activegenBefore, SADE.size() );
	    
     if(!this.Script.equals("NOFILE")){
		    
		    /*ADDING KNN FOR TEST FILE */
	        int trainRealClass[][];
	        int trainPrediction[][];
	            
	        trainRealClass = new int[datosTrain.length][1];
	        trainPrediction = new int[datosTrain.length][1];	
	       
	         nClases = SADE.getPosibleValuesOfOutput().size();
	
	       
	       //Working on training
	         int cont=0;
	        for (i=0; i<trainingDataSet.size(); i++) {
	             trainRealClass[i][0] = (int) trainingDataSet.get(i).getOutput(0);
	             trainPrediction[i][0] = evaluate(trainingDataSet.get(i).getInputs(),SADE.prototypeSetTodouble(), nClases, SADE.getClases(), 1);
	             
	             if(trainRealClass[i][0] ==  trainPrediction[i][0]){ cont++;}
	        }
	        
	        System.out.println("Acierto = "+ (cont*1.0)/(trainingDataSet.size()));
	        
	        Attribute entradas[];
	        Attribute salida;
	                     
	        entradas = Attributes.getInputAttributes();
	        salida = Attributes.getOutputAttribute(0);
	        String relation =  Attributes.getRelationName(); 
	            
	        writeOutput(this.ficheroSalida[0], trainRealClass, trainPrediction, entradas, salida, relation);
	        
	        int realClass[][] = new int[datosTest.length][1];
	        int prediction[][] = new int[datosTest.length][1];	
		
				
	        for (i=0; i<realClass.length; i++) {
		      realClass[i][0] = (int) testDataSet.get(i).getOutput(0);
		      prediction[i][0]= evaluate(testDataSet.get(i).getInputs(),SADE.prototypeSetTodouble(), nClases, SADE.getClases(), 1);
	        }
	            
	        writeOutput(this.ficheroSalida[1], realClass, prediction,  entradas, salida, relation);
	        saveRuntime(ficheroSalida[1]);
	    }
	    
  }
  
  
  /**
	 * Prints output files.
	 * 
	 * @param filename Name of output file
	 * @param realClass Real output of instances
	 * @param prediction Predicted output for instances
	 */
	public static void writeOutput(String filename, int [][] realClass, int [][] prediction, Attribute inputs[], Attribute output, String relation) {
	
		String text = "";
		
                
		/*Printing input attributes*/
		text += "@relation "+ relation +"\n";

		for (int i=0; i<inputs.length; i++) {
			
			text += "@attribute "+ inputs[i].getName()+" ";
			
		    if (inputs[i].getType() == Attribute.NOMINAL) {
		    	text += "{";
		        for (int j=0; j<inputs[i].getNominalValuesList().size(); j++) {
		        	text += (String)inputs[i].getNominalValuesList().elementAt(j);
		        	if (j < inputs[i].getNominalValuesList().size() -1) {
		        		text += ", ";
		        	}
		        }
		        text += "}\n";
		    } else {
		    	if (inputs[i].getType() == Attribute.INTEGER) {
		    		text += "integer";
		        } else {
		        	text += "real";
		        }
		        text += " ["+String.valueOf(inputs[i].getMinAttribute()) + ", " +  String.valueOf(inputs[i].getMaxAttribute())+"]\n";
		    }
		}

		/*Printing output attribute*/
		text += "@attribute "+ output.getName()+" ";

		if (output.getType() == Attribute.NOMINAL) {
			text += "{";
			
			for (int j=0; j<output.getNominalValuesList().size(); j++) {
				text += (String)output.getNominalValuesList().elementAt(j);
		        if (j < output.getNominalValuesList().size() -1) {
		        	text += ", ";
		        }
			}		
			text += "}\n";	    
		} else {
		    text += "integer ["+String.valueOf(output.getMinAttribute()) + ", " + String.valueOf(output.getMaxAttribute())+"]\n";
		}

		/*Printing data*/
		text += "@data\n";

		Files.writeFile(filename, text);
		
		if (output.getType() == Attribute.INTEGER) {
			
			text = "";
			
			for (int i=0; i<realClass.length; i++) {
			      
			      for (int j=0; j<realClass[0].length; j++){
			    	  text += "" + realClass[i][j] + " ";
			      }
			      for (int j=0; j<realClass[0].length; j++){
			    	  text += "" + prediction[i][j] + " ";
			      }
			      text += "\n";			      
			      if((i%10)==9){
			    	  Files.addToFile(filename, text);
			    	  text = "";
			      }     
			}			
			
			if((realClass.length%10)!=0){
				Files.addToFile(filename, text);
			}
		}
		else{
			
			text = "";
			
			for (int i=0; i<realClass.length; i++) {
			      
			      for (int j=0; j<realClass[0].length; j++){
			    	  text += "" + (String)output.getNominalValuesList().elementAt(realClass[i][j]) + " ";
			      }
			      for (int j=0; j<realClass[0].length; j++){
			    	  if(prediction[i][j]>-1){
			    		  text += "" + (String)output.getNominalValuesList().elementAt(prediction[i][j]) + " ";
			    	  }
			    	  else{
			    		  text += "" + "Unclassified" + " ";
			    	  }
			      }
			      text += "\n";
			      
			      if((i%10)==9){
			    	  Files.addToFile(filename, text);
			    	  text = "";
			      } 
			}			
			
			if((realClass.length%10)!=0){
				Files.addToFile(filename, text);
			}		
		}
		
	}//end-method
    
  
  /** 
	 * Calculates the Euclidean distance between two instances
	 * 
	 * @param instance1 First instance 
	 * @param instance2 Second instance
	 * @return The Euclidean distance
	 * 
	 */
	protected static double distance(double instance1[],double instance2[]){
		
		double length=0.0;

		for (int i=0; i<instance1.length; i++) {
			length += (instance1[i]-instance2[i])*(instance1[i]-instance2[i]);
		}
			
		length = Math.sqrt(length); 
				
		return length;
		
	} //end-method
    
	
	
	  /** 
	 * Calculates the Euclidean distance between two instances
	 * 
	 * @param instance1 First instance 
	 * @param instance2 Second instance
	 * @return The Euclidean distance
	 * 
	 */
	protected static double distanceWeighting(double instance1[],double instance2[], double Weights[]){
		
		double length=0.0;

		for (int i=0; i<instance1.length; i++) {
			length += ((instance1[i]-instance2[i])*(instance1[i]-instance2[i]))*Weights[i];
		}
			
		length = Math.sqrt(length); 
				
		return length;
		
	} //end-method
    
	
         
/** 
	 * Evaluates a instance to predict its class.
	 * 
	 * @param example Instance evaluated 
	 * @return Class predicted
	 * 
	 */
	public static int evaluate (double example[], double trainData[][],int nClasses,int trainOutput[],int k) {
	
		double minDist[];
		int nearestN[];
		int selectedClasses[];
		double dist;
		int prediction;
		int predictionValue;
		boolean stop;

		nearestN = new int[k];
		minDist = new double[k];
	
	    for (int i=0; i<k; i++) {
			nearestN[i] = 0;
			minDist[i] = Double.MAX_VALUE;
		}
		
	    //KNN Method starts here
	    
		for (int i=0; i<trainData.length; i++) {
		
		    dist = distance(trainData[i],example);

			if (dist > 0.0){ //leave-one-out
			
				//see if it's nearer than our previous selected neighbors
				stop=false;
				
				for(int j=0;j<k && !stop;j++){
				
					if (dist < minDist[j]) {
					    
						for (int l = k - 1; l >= j+1; l--) {
							minDist[l] = minDist[l - 1];
							nearestN[l] = nearestN[l - 1];
						}	
						
						minDist[j] = dist;
						nearestN[j] = i;
						stop=true;
					}
				}
			}
		}
		
		//we have check all the instances... see what is the most present class
		selectedClasses= new int[nClasses];
	
		for (int i=0; i<nClasses; i++) {
			selectedClasses[i] = 0;
		}	
		
		for (int i=0; i<k; i++) {
             //      System.out.println("nearestN i ="+i + " =>"+nearestN[i]);
              // System.out.println("trainOutput ="+trainOutput[nearestN[i]]);
                
			selectedClasses[trainOutput[nearestN[i]]]+=1;
		}
		
		prediction=0;
		predictionValue=selectedClasses[0];
		
		for (int i=1; i<nClasses; i++) {
		    if (predictionValue < selectedClasses[i]) {
		        predictionValue = selectedClasses[i];
		        prediction = i;
		    }
		}
		
		return prediction;
	
	} //end-method	
	

	

  public void leerConfiguracion (String ficheroScript) {

    String fichero, linea, token;
    StringTokenizer lineasFichero, tokens;
    byte line[];
    int i, j;

    ficheroSalida = new String[2];

    if(ficheroScript.equals("NOFILE")){
    	System.out.println("There is no configuration file: Applying Auto-parameters");
    	    	
    	ficheroSalida[0] = "salida.dat";
    	ficheroSalida[1] = "otro.dat";
    	ficheroTraining = "intermediate.dat";
    	tamPoblacion = 30;
    	nEvalSSMA = 10000;
    	pCross = 0.5;
    	pMut = 0.001;
    	kNeigh = 1;
    	distanceEu = true;
    	PopulationSize = 40;
    	this.maxEval_IG = 100;
	    this.iterSFGSS = 8;
	    this.iterSFHC = 20;
	    this.Fl =  0.1;
	    this.Fu =  0.9;
	     tau = new double[4];
	    this.tau[0] =  0.1;
	    this.tau[1] =  0.1;
	    this.tau[2] =  0.03;
	    this.tau[3] =  0.07;;
	    this.Strategy = 3;

    	
    }else{
	    fichero = Fichero.leeFichero (ficheroScript);
	    
	
	    lineasFichero = new StringTokenizer (fichero,"\n\r");
	
	    lineasFichero.nextToken();
	    linea = lineasFichero.nextToken();
	
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    token = tokens.nextToken();
	
	    /*Getting the name of training and test files*/
	    line = token.getBytes();
	    for (i=0; line[i]!='\"'; i++);
	    i++;
	    for (j=i; line[j]!='\"'; j++);
	    ficheroTraining = new String (line,i,j-i);
	    
		for (i=j+1; line[i]!='\"'; i++);
		i++;
		for (j=i; line[j]!='\"'; j++);
		ficheroValidation = new String (line,i,j-i);
	    
	    for (i=j+1; line[i]!='\"'; i++);
	    i++;
	    for (j=i; line[j]!='\"'; j++);
	    ficheroTest = new String (line,i,j-i);
	
	    
	    //Parameters.assertBasicArgs(ficheroSalida);
	    
	    
	
	    
	    /*Obtainin the path and the base name of the results files*/
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    token = tokens.nextToken();
	
	    /*Getting the name of output files*/
	    line = token.getBytes();
	    for (i=0; line[i]!='\"'; i++);
	    i++;
	    for (j=i; line[j]!='\"'; j++);
	    ficheroSalida[0] = new String (line,i,j-i);
	    for (i=j+1; line[i]!='\"'; i++);
	    i++;
	    for (j=i; line[j]!='\"'; j++);
	    ficheroSalida[1] = new String (line,i,j-i);
	
	    /*Getting the seed*/
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    semilla = Long.parseLong(tokens.nextToken().substring(1));
	
	    /*Getting the size of the poblation and the number of evaluations*/
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    tamPoblacion = Integer.parseInt(tokens.nextToken().substring(1));
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
        
        // ISAAC: we need to modify this, to multiply nEval * numberFeatures * number of instances.
	    nEvalSSMA = Double.parseDouble(tokens.nextToken().substring(1));
	
	    /*Getting the probabilities of evolutionary operators*/
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    pCross = Double.parseDouble(tokens.nextToken().substring(1));
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    pMut = Double.parseDouble(tokens.nextToken().substring(1));
	
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    kNeigh = Integer.parseInt(tokens.nextToken().substring(1));
	 
	    /*Getting the type of distance function*/
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    distanceEu = tokens.nextToken().substring(1).equalsIgnoreCase("Euclidean")?true:false;    
	    
	    
	    
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    this.PopulationSize = Integer.parseInt(tokens.nextToken().substring(1));
	    
        // isaac suggestion:  to remove maxEval_IG.
	    //As this variable will be modify later in code, we just read the val from the config file 
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    this.maxEval_IG = Integer.parseInt(tokens.nextToken().substring(1));
	    
	    
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    this.iterSFGSS = Integer.parseInt(tokens.nextToken().substring(1));
	    
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    this.iterSFHC = Integer.parseInt(tokens.nextToken().substring(1));
	    
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    this.Fl =  Double.parseDouble(tokens.nextToken().substring(1));
	    
	    
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    this.Fu =  Double.parseDouble(tokens.nextToken().substring(1));
	    
	    tau = new double[4];
	    
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    this.tau[0] =  Double.parseDouble(tokens.nextToken().substring(1));
	   
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    this.tau[1] =  Double.parseDouble(tokens.nextToken().substring(1));
	    
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    this.tau[2] =  Double.parseDouble(tokens.nextToken().substring(1));
	    
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    this.tau[3] =  Double.parseDouble(tokens.nextToken().substring(1));
	    
	    
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    this.Strategy = Integer.parseInt(tokens.nextToken().substring(1));
	    
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    this.perc = Integer.parseInt(tokens.nextToken().substring(1));
	    
	    
	 //   System.out.print("\nIsaac dice:  tau3"+this.tau[3] +"\n");
	  
    
    }
    
  }
}
