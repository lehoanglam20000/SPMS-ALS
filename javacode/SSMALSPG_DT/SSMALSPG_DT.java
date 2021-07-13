//
//  SSMA.java
//
//  Salvador Garc�a L�pez
//
//  Created by Salvador Garc�a L�pez 3-10-2005.
//  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
//

package keel.Algorithms.Instance_Generation.SSMALSPG_DT;

import keel.Algorithms.Preprocess.Basic.*;

import keel.Algorithms.Instance_Generation.Basic.Prototype;
import keel.Algorithms.Instance_Generation.Basic.PrototypeGenerationAlgorithm;
import keel.Algorithms.Instance_Generation.Basic.PrototypeSet;
import keel.Algorithms.Instance_Generation.Basic.PublicVar;
import keel.Algorithms.Instance_Generation.Basic.UtilFunc;
import keel.Algorithms.Instance_Generation.SSMALSPG.Cromosoma;
import keel.Algorithms.Instance_Generation.utilities.*;
import keel.Algorithms.Instance_Generation.Basic.PrototypeGenerator;

import java.util.*;

import org.core.*;

import umontreal.ssj.functionfit.LeastSquares;

import keel.Dataset.Attribute;
import keel.Dataset.Attributes;
import keel.Dataset.InstanceAttributes;
import keel.Dataset.InstanceSet;
import umontreal.ssj.functionfit.LeastSquares;

import org.core.*;

import java.util.StringTokenizer;
import java.util.ArrayList;
import java.util.Arrays;

public class SSMALSPG_DT extends Metodo {

  /*Own parameters of the algorithm*/

  private long semilla;
  private int tamPoblacion;
  private double nEvalSSMA;
  private double pCross;
  private double pMut;
  private int kNeigh;
  public String Script; // para releer par�metros..
  private PrototypeSet trainingDataSet;
  private PrototypeSet testDataSet;
  private PrototypeGenerator generador;
  
  //Par�metros for LSPG
  private int k;
  private int perc; 
  private double rho;
  private String mode; // individual creation mode:  cap, toroidal or reflection.
  private int maxEval_IG;
  private int MaxEvalPerDimension;
  private boolean Reflection = false; // variable to control if reflection needs to happen.
  
  protected int numberOfClass;


  protected int numberOfPrototypes;  // Particle size is the percentage
  protected int numberOfStrategies; // number of strategies in the pool
  protected double [][] distMatrix = null;
  protected double [] backup_col = null;
  private boolean suffleIndex;  //The mode will decide how to create a SA: select the instances randomly or iteratively
  public SSMALSPG_DT (String ficheroScript) {
	super (ficheroScript);
	suffleIndex = PublicVar.ISRANDOM;
  }


  public SSMALSPG_DT(String ficheroScript, InstanceSet train) {
		super (ficheroScript, train);
		suffleIndex = PublicVar.ISRANDOM;
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
  
      double acc = 100.0* (wellClassificated / (double)test.size());
      return acc;
  }
  
  
  
  public double[] applyLeastSquare(Pair<ArrayList<PrototypeSet>,ArrayList<Double>> Surr){
	  	  
	  double[][] points = new double[Surr.first().size()][];
	  double[] fitnesses = new double[Surr.first().size()];

	  ArrayList<Double> helper =  (ArrayList<Double>) Surr.second();
	    
	  Surr.first().get(0).printSet();
	  
	  System.out.println("Points : "+ Surr.first().size());
	  
	  for(int i=0; i<Surr.first().size(); i++){ 
		  points[i] = Surr.first().get(i).prototypeSetTodoubleVector();
		//  System.out.println("points[i].length: "+points[i].length);
		
		  fitnesses[i] = helper.get(i).doubleValue();
		  //System.out.println(fitnesses[i]);
	  }
	  
	  
	 // (Y.length <= X[0].length + 1

	  System.out.println(fitnesses.length);
	  System.out.println(points[0].length+1);

	  
	 return LeastSquares.calcCoefficients(points, fitnesses);
  }
  
  public double FitnessApproximation(double[] coefficients, PrototypeSet sol){
	  double fitness = 0;
	  for (int i = 0; i<coefficients.length; i++){
		  int proto = i/sol.size();
		  int feature= i%sol.size();
		  
		  fitness+= coefficients[i]*sol.get(proto).getInput(feature);
	  }
	  return fitness;
  }
  
  /**
   * Functon to compute the fitness using corrections for nominal values.
   * @param sol
   * @return
   */
  
  public double fitnessComputation(PrototypeSet sol){
	  PrototypeSet nominalPopulation = new PrototypeSet();
      nominalPopulation.formatear(sol);
	  return generador.accuracy(nominalPopulation,trainingDataSet);
  }
  
  
  public double CheckBoundaries(double newValue, double currentValue){
	  if(this.mode.equalsIgnoreCase("cap")){
		  if (newValue>1.0) newValue = 1.0;
		  else if(newValue<0.0) newValue =0.0;
	  }else if(this.mode.equalsIgnoreCase("toroidal")){
		  if (newValue>1.0) newValue = newValue-1.0;
		  else if(newValue<0.0){ 
			  newValue =1.0+newValue;
		  //	  System.out.println("here I am: "+newValue);
		  }
	  }else{
		 // System.out.println("\n reflection \n:"+this.mode+";");
		  if(newValue<0.0){
			  newValue = currentValue+this.rho;
		  }else if(newValue>1.0) newValue = currentValue -this.rho/2;
		  
		  this.Reflection = true;

	  }
	  return newValue;
  }
  
  public void buildDistanceMatrix(PrototypeSet sol){
	  
	  if (backup_col == null)
		  backup_col = new double[trainingDataSet.size()];
	  if (distMatrix == null) {
		  distMatrix = new double[trainingDataSet.size()][];
		  for (int i = 0; i<trainingDataSet.size(); i++) {
			  distMatrix[i] = new double[sol.size()];
		  }
	  }
		 
	  PrototypeSet nominalPopulation = new PrototypeSet();
      nominalPopulation.formatear(sol);
      int p_index = 0;
      for(Prototype p_InTR : trainingDataSet)
      {
          for (int i=0; i<nominalPopulation.size(); i++)
          {
              Prototype prottp_InRS = nominalPopulation.get(i);
              double currDist = Distance.euclideanDistance(p_InTR, prottp_InRS);
              if (currDist > 0)
            	  distMatrix[p_index][i] = currDist;
              else
            	  distMatrix[p_index][i] = Double.MAX_VALUE;
          }
          p_index++;
      }
  } 
  /*
   * idx is the index of the inst in sol that we are modifying its features.
   * if we know what prototype we are editing. for example: idx i
	  we update the distance matrix, then get the label of the insntace that has the smallest distance
	  
	  
	   
   */
  public double fitnessDelta (PrototypeSet sol, int idx){
	  
	  PrototypeSet nominalPopulation = new PrototypeSet();
      nominalPopulation.formatear(sol);
      if (idx == PublicVar.NONMODIFIED) {
		  return computeAcc(nominalPopulation); 
	  }
      else {
		  Prototype modifidedInst = nominalPopulation.get(idx);
	      double currDist;
	      
	      for (int i=0; i<trainingDataSet.size(); i++)
	      {
	          Prototype prottp = trainingDataSet.get(i);
	          currDist = Distance.euclideanDistance(prottp, modifidedInst);
	          if (currDist == 0) {
	        	 // System.out.println("Use leave one out validation");
	          }
	          if (currDist > 0) {
	        	  distMatrix[i][idx] = currDist;
	          }
	          else //currDist ==0: (cannot have currDist <0): We apply leave-one-out validation, so not select this one 
	        	  //as the nearest neighbour for prediction an unseen sample. 
	        	  distMatrix[i][idx] = Double.MAX_VALUE;
	           
	      }
			      
		  double acc = computeAcc(nominalPopulation);
		  return acc;
      }
	  
	  /*
	 double acc2 = fitnessComputation_check(sol);
	 if (acc != acc2) {
		 System.err.println("ERRor of computing fitness");
	 }
	 System.err.println("Acc1 - acc2 " + acc + " -- " + acc2);
	  */
      
  }
  
 
  public double computeAcc(PrototypeSet nominalPopulation) {
	  double wellclassified = 0;
	  for(int i = 0; i<distMatrix.length; i++) {
		  Prototype p = trainingDataSet.get(i);
		  int indexNN = UtilFunc.getIndexOfMinDist(distMatrix[i]);
		  Prototype nearestNeighbor = nominalPopulation.get(indexNN);
          
          if(p.getOutput(0) == nearestNeighbor.getOutput(0)) {
        	  wellclassified++;
          }
         
	  }
	  double acc = 0;
	  
	  acc = wellclassified / distMatrix.length * 100;
	  return acc;
  }
  /**
   * Generate a reduced prototype set by the LSPG generator method.
   */
  

  public PrototypeSet reduceSet(PrototypeSet initial)
  {
  System.out.print("\nThe algorithm  LSPG is starting...\nComputing...\n");
	  numberOfPrototypes = initial.size();
	  System.out.println("Number of prototypes, result set = "+numberOfPrototypes+ "\n");
	  if(numberOfPrototypes < trainingDataSet.getPosibleValuesOfOutput().size()){
		  System.out.println("Number of prototypes less than the number of clases");
		  numberOfPrototypes = trainingDataSet.getPosibleValuesOfOutput().size();
	  }
	  generador = new PrototypeGenerator(trainingDataSet);
	  this.numberOfClass = trainingDataSet.getPosibleValuesOfOutput().size();
	  /*
	   * Parameters:
	   * 	perc: [0,100] - percentage of reduction from original set: 
	   *    rho: [0,0.4] - constant value that indicates the movement on an attribute; it assumed the solution attributes are in [0,1].
	   *    mode: {cap, toroidal,reflection} - when creating a perturbation and goes over the range.. cap will simple not allow greater or less than 0,1; toroidal will go around the range  and relfection will change the sign of the perturbations.
	   *    MaxEval should be set to 5000*N..
	   * */
	  
	  // Step 0: random initialisation of a potential selection of prototypes from TR; 
	  // size of the solution depends on the input parameter `numberOfPrototypes`.
	  
	  PrototypeSet solution = initial.clone();
	  int dimensions = initial.get(0).numberOfInputs();
	  
	  System.out.println("dimensions: "+dimensions + ", MaxEvalPerDimension"+ MaxEvalPerDimension);
	 //Psxhll note: As SSMA is used in a combination with other instance generation techniques. The number of evaluation is modified.
	  this.maxEval_IG = (int) (trainingDataSet.size() * (dimensions*PublicVar.EvalWeight - PublicVar.SSMAWeight));
      System.out.print("\nIsaac checks:  " + k + " mode= "+mode+ " rho= "+ rho+"; max eval = "+this.maxEval_IG+"\n");

	  
	  // Ensure that we have at least a representative of each class.
	  PrototypeSet clases[] = new PrototypeSet [this.numberOfClass];
	  for(int i=0; i< this.numberOfClass; i++){
		  clases[i] = new PrototypeSet(trainingDataSet.getFromClass(i));
	  }
	
	  for(int j=0; j< this.numberOfClass; j++){
		  if(solution.getFromClass(j).size() ==0 && clases[j].size()!=0){
			  solution.add(clases[j].getRandom());
		  }
	  }
	  
	  int numFeatures = solution.get(0).numberOfInputs();
	  // main while loop.
	  int totalEval = 0;
	  int sameFitness=0;
	  int instanceIdx =-1 ;
	  buildDistanceMatrix(solution);
	  while (totalEval < this.maxEval_IG & this.rho > 0.0001){
		  //PrototypeSet OriginalSolution = solution.clone(); // should this be here or within while loop?

		  int numEvals = 0;
		  
		  double fitness = fitnessDelta(solution, PublicVar.NONMODIFIED);   // no need to format as nomimal Population as they are originally selected from TR.
		  
		  if(totalEval==0)
			  System.out.println(totalEval+"\t"+fitness);

		  if(fitness == 100){
			  System.out.println("\n ****** \n optimal value for the fitness");
			  break;
		  }
		  
		  numEvals++;
		  // solution.printSet();
		  // System.out.println("Fitness of the initial solution: "+fitness+ "\n");
		  		 
		  // we are going to treat the solution as a 1D - vector.
		  int sizeVector = solution.size()*solution.get(0).numberOfInputs();
		 // System.out.println("Size of the vector: "+sizeVector+ "\n");
		  
		  // version 1:  iterative... from instance 0 to instance N
		  // version 2: randomly selected instances. 
		  
		  //ArrayList<Integer> indexes =  RandomGenerator.generateDifferentRandomIntegers(0, sizeVector-1);
		  
		 // for (int i=0; i<sizeVector; i++){
		  int randomIndexes [] = UtilFunc.createShuffledIndexArr(sizeVector);
		  int goodMoves=0;
		//  System.out.println("\n Iteration: \n");
		  double currentFitness = fitness;
		  for (int iter=0;iter<sizeVector; iter++){	  
			  this.Reflection = false;

			  // determine the position to modify depending on mode:
			  //int pos =  iter;
			  int pos = (suffleIndex == false) ? iter : randomIndexes[iter]; //Select the position using a random or iterative
			  
			  // calculate prototype and feature.
			  int proto = pos/numFeatures;
			  int feature= pos%numFeatures;
			  instanceIdx = pos / solution.get(0).numberOfInputs();
			//  System.out.println("pos:  "+pos+ "Proto "+proto + ", feature: "+feature);
			  double currentValue = solution.get(proto).getInput(feature);
			  //double newValue = (RandomGenerator.Rand()<0.5) ? currentValue+this.rho : currentValue-this.rho;
			  double newValue =  currentValue-this.rho;
			
			//back up the col before modifying
		      for (int i=0; i<trainingDataSet.size(); i++)
		      {
		    	  backup_col[i] = distMatrix[i][instanceIdx];
		      }
			  // Check boundaries:
			  newValue= CheckBoundaries(newValue, currentValue);
			  solution.get(proto).setInput(feature, newValue);
			  // solution.get(proto).applyThresholds(); // check we didn't go above or below boundaries.
			  //Check fitness of this move
			  double newFitness = fitnessDelta(solution, instanceIdx); 
			  numEvals++;
			  
			  if(newFitness>=fitness){
				  // System.out.print("Good move on " + iter + " ;" + newValue+"; realNV : "+ realNV+" \t");
				  if(newFitness>fitness){
					  int currentEval = totalEval+numEvals;
					  System.out.println(currentEval+"\t"+newFitness);
				  }
				  fitness = newFitness; // if that was a good move, we update the fitness
				  goodMoves++;
			  }else{
				  // we are going to try another correction:
				  if(!Reflection)
					  newValue = currentValue + this.rho/2.0;
				  else
					  newValue = currentValue - this.rho/2.0;
				  
				  // Check boundaries:
				//  realNV = newValue;
				  newValue= CheckBoundaries(newValue, currentValue);
				  
				  solution.get(proto).setInput(feature, newValue);
				  
				 // solution.get(proto).applyThresholds(); // check we didn't go above or below boundaries.
				  
				  
				  newFitness = fitnessDelta(solution, instanceIdx); 
				  numEvals++;
				  		  
				  if(newFitness>=fitness){ // if we have decreased the fitness... we restore it
					  if(newFitness>fitness){
						  int currentEval = totalEval+numEvals;
						  System.out.println(currentEval+"\t"+newFitness);
					  }
					  fitness = newFitness;
					  goodMoves++;
					//  System.out.print("Good move on (half rho+) " + iter + " ;" + newValue+"; realNV : "+ realNV+" \t");

				  }else{
					  solution.get(proto).setInput(feature,currentValue); // restore value
					//need to update the distMatrix, reverse to the state when the feature was not changed. 
				      for (int i=0; i<trainingDataSet.size(); i++)
				      {
				    	  distMatrix[i][instanceIdx] =backup_col[i];
				      }
				  }
	
			  }
			  
			  //System.out.println("fitness: "+fitness);

			 
		  }
		  
		  if(currentFitness == fitness){
			  sameFitness++;
		  }
		//  System.out.println("\nI've made "+goodMoves+ " good moves; same fitness: "+sameFitness);// in " + sizeVector+" options; fitness: "+ fitness );
	
	
		  // Step 5: check if we have done something with solution.
		  //if (solution.equals(OriginalSolution)){ //TODO: this needs testing
		  if (goodMoves ==0 | sameFitness > 10){
			  sameFitness =0;
		  	  this.rho/=2;
		  	/*
			  	if (goodMoves == 0)
			  		System.out.println("I didn't change anything , reducing rho: "+this.rho);
			  	else
			  		System.out.println("I was in plateau, reducing rho "+this.rho );
			  		*/
		  }
		  
		  totalEval +=  numEvals;
;
		  
	  }// end  main-while loop
	  
	  	int evalSaved = (int)(this.maxEval_IG-totalEval>0? this.maxEval_IG-totalEval: 0);
	   fitCollect.collectEvalSaved(this.maxEval_IG, evalSaved); //Because while loop can stop with rho conditions, 
	   //We may save a number of (this.maxEval_IG-totalEval) evaluations.
        System.out.println("Finish - LSPG");
        
        System.out.println("I should have done a max evaluations: "+this.maxEval_IG);
        System.out.println("Total number of evaluations: "+totalEval);
		PrototypeSet nominalPopulation = new PrototypeSet();
        nominalPopulation.formatear(solution);
		System.err.println("\n% de acierto en training Nominal " + classficationAccuracy1NN(nominalPopulation,trainingDataSet));
				
		return nominalPopulation;
  }
  
  
  /* MEzcla de algoritmos */
  public void ejecutar () {

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

    /*Initial evaluation of the population*/
    for (i=0; i<tamPoblacion; i++) {
      poblacion[i].evaluacionCompleta(nClases, kNeigh, clasesTrain);
    }

    umbralOpt = 0;

    /*Until stop condition*/
    //Psxhll note: As SSMA is used in a combination with other instance generation techniques. The number of evaluation is modified.
    if (PublicVar.OriginalSetting) {
    	nEvalSSMA = PublicVar.OriginalEvalSSMA;
    }else {
    	nEvalSSMA = PublicVar.SSMAWeight * datosTrain.length;
    }
    while (ev < nEvalSSMA) {
	    	Arrays.sort(poblacion);
	    	//System.out.println("Ev = " + ev);
		         
	
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
	      hijos[0].evaluacionCompleta(nClases, kNeigh, clasesTrain);
	      hijos[1].evaluacionCompleta(nClases, kNeigh, clasesTrain);
	      ev+=2;
	      GAeffort += 2;
	      temporal = ev;
	      
	      //Psxhll: If newly found sol is better than the worst solution in the pop, perform a local search
	      if (hijos[0].getFitness() > poblacion[tamPoblacion-1].getFitness() || Randomize.Rand() < 0.0625) {
	    	  ev += hijos[0].optimizacionLocal(nClases, kNeigh, clasesTrain,dMatrix,umbralOpt, datosTrain, realTrain, nominalTrain, nulosTrain, distanceEu);
	      }
	      
	      if (hijos[1].getFitness() > poblacion[tamPoblacion-1].getFitness() || Randomize.Rand() < 0.0625) {
	          ev += hijos[1].optimizacionLocal(nClases, kNeigh, clasesTrain,dMatrix,umbralOpt, datosTrain, realTrain, nominalTrain, nulosTrain, distanceEu);
	      }
	
	      LSeffort += (ev - temporal);
	
	      /*Replace the two worst*/
	      if (hijos[0].getFitness() > poblacion[tamPoblacion-1].getFitness()) {
	        poblacion[tamPoblacion-1] = new Cromosoma (kNeigh, datosTrain.length, hijos[0]);
	      }
	      if (hijos[1].getFitness() > poblacion[tamPoblacion-2].getFitness()) {
	        poblacion[tamPoblacion-2] = new Cromosoma (kNeigh, datosTrain.length, hijos[1]);
	      }
	
	/*      System.out.println(ev + " - (" + umbralOpt + ")" +
	                         (double)poblacion[0].getFitnessAc()*100.0/(double)datosTrain.length + " / " +
	                         (double)poblacion[tamPoblacion-1].getFitnessAc()*100.0/(double)datosTrain.length + " - " +
	                         (1.0-((double)poblacion[0].genesActivos()/(double)datosTrain.length))*100.0 + " / " +
	                         (1.0-((double)poblacion[tamPoblacion-1].genesActivos()/(double)datosTrain.length))*100.0);*/
	      
		}//End while loop
    
	    
	    
		Arrays.sort(poblacion);
		nSel = poblacion[0].genesActivos();
	    //System.out.println("\nKNN accuracy from the reduced set of SSMA "+ poblacion[0].getFitness());
	    //poblacion[0].evaluacionCompleta(nClases, kNeigh, clasesTrain);
	    //System.out.println("\nKNN accuracy from the reduced set of SSMA "+ poblacion[0].getFitness());
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
	    System.out.println("--- Runtime of SSMA  for "+ relation + " " + (double)(System.currentTimeMillis()-tiempo)/1000.0 + "s");
	
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
	   // trainingDataSet.print();
	   //this.numberOfPrototypes = (int)Math.floor((trainingDataSet.size())*ParticleSize/100.0);
	   // training.print();
		 
	
	    //--------------------------------------------------------
			// If SSMA fails (It may happen, I know!):
   		 	  
	  if(training.size() <2){
		  System.out.println("This fold has training.size() < 2, create a random solution");
		  this.numberOfPrototypes = (int)Math.round(trainingDataSet.size() * this.perc/100);
		  generador = new PrototypeGenerator(trainingDataSet);
		  training = generador.selecRandomSet(numberOfPrototypes, true).clone() ;
		  // red .95
		  
		  // Aseguro que al menos hay un representante de cada clase.
		  this.numberOfClass = trainingDataSet.getPosibleValuesOfOutput().size();
		  PrototypeSet clases[] = new PrototypeSet [this.numberOfClass];
		  for(i=0; i< this.numberOfClass; i++){
			  clases[i] = new PrototypeSet(trainingDataSet.getFromClass(i));
			  
			 // System.out.println("Clase "+i+", size= "+ clases[i].size());
		  }
		
		  for(i=0; i< training.size(); i++){
			  for(j=0; j< this.numberOfClass; j++){
				  if(training.getFromClass(j).size() ==0 && clases[j].size()!=0){
					  
					  training.add(clases[j].getRandom());
				  }
			  }
		  }
	  }
	  //--------------------------------------------------------
     PrototypeSet SADE = reduceSet(training); 
     //Psxhll: Print out the fitness collected in child class. The child class needs to override the method in the parent class
	    fitCollect.writeToFile(ficheroSalida[0]);
	    fitCollect.writeToFile_ReductionRate(ficheroSalida[0], training.size(), trainingDataSet.size());
     //SADE.save(ficheroSalida[0]); 
     //    SADE.print();

     //Copy the test input file to the output test file
     // KeelFile.copy(inputFilesPath.get(TEST), outputFilesPath.get(TEST));
       
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
	        
	        System.out.println("Time elapsed for the entire SSMALSPG:" + (double)(System.currentTimeMillis()-tiempo)/1000.0 + "s");
	           
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
    	mode="toroidal";
    	rho=0.4;
    	MaxEvalPerDimension=100;
    	

    	
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
	    nEvalSSMA = Double.parseDouble(tokens.nextToken().substring(1));
	    //Psxhll note: As SSMA is used in a combination with other instance generation techniques. The number of evaluation is modified. 
	    
	
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
	    
	    
	     // this.MaxEval = parameters.getNextAsInt()*trainingDataSet.get(0).numberOfInputs()*numberOfPrototypes;
	      
	      
	    
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    this.mode = tokens.nextToken().substring(1);
	    
	    
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    this.rho = Double.parseDouble(tokens.nextToken().substring(1));
	    
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    this.MaxEvalPerDimension = Integer.parseInt(tokens.nextToken().substring(1));
	    
	    linea = lineasFichero.nextToken();
	    tokens = new StringTokenizer (linea, "=");
	    tokens.nextToken();
	    this.perc = Integer.parseInt(tokens.nextToken().substring(1));
	
	    
	 //   System.out.print("\nIsaac dice:  tau3"+this.tau[3] +"\n");

    
    }
    
  }
}

