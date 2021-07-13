
/*
	LSPG.java
	Isaac Triguero Velazquez and Ferrante Neri
	
	October 2019

*/

package keel.Algorithms.Instance_Generation.LSPG_DT;

import keel.Algorithms.Instance_Generation.Basic.PrototypeSet;
import keel.Algorithms.Instance_Generation.Basic.PublicVar;
import keel.Algorithms.Instance_Generation.Basic.UtilFunc;
import keel.Algorithms.Instance_Generation.Basic.PrototypeGenerator;
import keel.Algorithms.Instance_Generation.Basic.Prototype;
import keel.Algorithms.Instance_Generation.Basic.PrototypeGenerationAlgorithm;
import java.util.*;

import keel.Algorithms.Instance_Generation.utilities.*;
import keel.Algorithms.Instance_Generation.utilities.KNN.*;

import org.core.*;

import umontreal.ssj.functionfit.LeastSquares;


/**

 * @author Isaac Triguero and Ferrante Neri
 * @version 1.0
 */
public class LSPG_DTGenerator extends PrototypeGenerator {

  /*Own parameters of the algorithm*/
  
  // We need the variable K to use with k-NN rule
  private int k;
 
 
  private double rho;
  private String mode; // individual creation mode:  cap, toroidal or reflection.
  private int MaxEval;
  private boolean suffleIndex;  //The mode will decide how to create a SA: select the instances randomly or iteratively
  private boolean Reflection = false; // variable to control if reflection needs to happen.
  protected int numberOfClass;
  protected int numberOfPrototypes;  // individual size is the percentage
  protected double [][] distMatrix = null;
  protected double [] backup_col = null;

  
  
  /**
   * Build a new LSPGGenerator Algorithm
   * @param t Original prototype set to be reduced.
   * @param perc Reduction percentage of the prototype set.
   */
  
  public LSPG_DTGenerator(PrototypeSet _trainingDataSet, int neigbors, int perc, String searchMode, double rhoStep, int eval)
  {
      super(_trainingDataSet);
      algorithmName="LSPG";
      
      this.k = neigbors;
      this.mode = searchMode;
      this.rho = rhoStep;
      this.MaxEval = eval;
      this.numberOfPrototypes = getSetSizeFromPercentage(perc);
      suffleIndex = PublicVar.ISRANDOM;
  }
  
  

  /**
   * Build a new LSPGGenerator Algorithm
   * @param t Original prototype set to be reduced.
   * @param params Parameters of the algorithm (only % of reduced set).
   */
  public LSPG_DTGenerator(PrototypeSet t, Parameters parameters)
  {
      super(t, parameters);
      algorithmName="LSPG";
      this.k =  parameters.getNextAsInt();
      this.mode = parameters.getNextAsString().trim();
      this.rho = parameters.getNextAsDouble();      
      int perc =  parameters.getNextAsInt();
      this.numberOfPrototypes = getSetSizeFromPercentage(perc);

      //this.MaxEval = parameters.getNextAsInt()*trainingDataSet.get(0).numberOfInputs()*numberOfPrototypes;
      this.MaxEval = (int)(PublicVar.EvalWeight*trainingDataSet.get(0).numberOfInputs()*trainingDataSet.size());
      suffleIndex = PublicVar.ISRANDOM;
      this.numberOfClass = trainingDataSet.getPosibleValuesOfOutput().size();
      
      System.out.print("\nIsaac checks:  " + k + " mode= "+mode+ " perc=  "+ perc + " rho= "+ rho+"; max eval = "+this.MaxEval+"\n");
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
  
  /*
   *Option 2:  Use the distance matrix, and another array to store the min distance and index of mindist

  */
  public double fitnessDelta_option2 (PrototypeSet sol, int idx, double fitness){
 	  
 	  PrototypeSet nominalPopulation = new PrototypeSet();
      nominalPopulation.formatear(sol);
      if (idx == PublicVar.NONMODIFIED) {
 		  return computeAcc(nominalPopulation); 
 	  }
      else {
    	  double backup[][] = UtilFunc.copyArr2D(distMatrix);
 		  Prototype modifidedInst = nominalPopulation.get(idx);
 	      double currDist;
 	      
 	      for (int i=0; i<trainingDataSet.size(); i++)
 	      {
 	          Prototype prottp = trainingDataSet.get(i);
 	          currDist = Distance.euclideanDistance(prottp, modifidedInst);
 	          
 	          
        	  if (currDist > 0) {
	        	  distMatrix[i][idx] = currDist;
	          }
	          else //currDist ==0: (cannot have currDist <0): We apply leave-one-out validation, so not select this one 
	        	  //as the nearest neighbour for prediction an unseen sample. 
	        	  distMatrix[i][idx] = Double.MAX_VALUE;
 	          
 	         // System.out.println(backup[0][0] + "  -  " + backup[0][1]);
 	      }
 	      //Before changing, label is decided by instance 1, now we modify instance 1, the distance is updated to the largestWhen currDist > backup[i][0], still change, now the label of unseen sample is determined by the second samples
 		  double acc = computeAcc(nominalPopulation);
 		  if (acc >= fitness) {
 			  distMatrix = UtilFunc.copyArr2D(backup); 
 		  }
 		
 		  return acc;
      }
  }
 
 
  public double computeAcc(PrototypeSet nominalPopulation) {
	 // UtilFunc.print2DArr(distMatrix);
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
	  //System.out.printf("Fitness of this matrix: %.4f\n ", acc);
	  return acc;
  }
  public double fitnessComputation_check(PrototypeSet sol){
	  PrototypeSet nominalPopulation = new PrototypeSet();
      nominalPopulation.formatear(sol);
	  return accuracy(nominalPopulation,trainingDataSet);
  }
  public double fitnessComputation(PrototypeSet sol){
	 
	  PrototypeSet nominalPopulation = new PrototypeSet();
      nominalPopulation.formatear(sol);
      
      int wellClassificated = 0;
      for(Prototype p : trainingDataSet)
      {
          int indexNN = _1nn(p, nominalPopulation);
          Prototype nearestNeighbor = nominalPopulation.get(indexNN);
          if(p.getOutput(0) == nearestNeighbor.getOutput(0))
              ++wellClassificated;
      }
      double acc = 100.0 * (wellClassificated / (double)trainingDataSet.size());
      
      return acc;
      
  }
 	
  
  
  public int _1nn(Prototype inst, PrototypeSet dataSet)
  {
      int indexNN = 0;
      double minDist =Double.POSITIVE_INFINITY;
      double currDist;
      for (int i=0; i<dataSet.size(); i++)
      {
          Prototype prottp = dataSet.get(i);
          currDist = Distance.euclideanDistance(prottp, inst);
           if(currDist >0){
              if (currDist < minDist)
              {
                  minDist = currDist;
                  indexNN =i;
              }
          }
      }
      
      return indexNN;
  }
  
  public static int[] _2nn(Prototype current, PrototypeSet dataSet)
  {
	  int k = 2;
      boolean stop;
      int nearestN [] = new int[k];
      double minDist[] = new double[k];
	  double dist;
	  for (int i=0; i<k; i++) {
		nearestN[i] = 0;
		minDist[i] = Double.MAX_VALUE;
	  }
      for (int i=0; i<dataSet.size(); i++)
      {
          Prototype pi = dataSet.get(i);
          dist = Distance.euclideanDistance(pi,current);
           if(dist >0){
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
      return nearestN;
  }
  protected static double euclideanDistance(double instance1[],double instance2[]){
		
		double length=0.0;

		for (int i=0; i<instance1.length; i++) {
			length += (instance1[i]-instance2[i])*(instance1[i]-instance2[i]);
		}
			
		length = Math.sqrt(length); 
				
		return length;
		
	} //end-method
	
  
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
  
  public PrototypeSet reduceSet_OriginalSetting() {
	  this.MaxEval = PublicVar.OriginalEvalLSPG *trainingDataSet.get(0).numberOfInputs()*numberOfPrototypes; //100 * noOFFeautres on the solution. 
      return reduceSet();
  }
  /**
   * Generate a reduced prototype set by the LSPGGenerator method.
   * @return Reduced set by LSPGGenerator's method.
   */
    
  public PrototypeSet reduceSet()
  {
	  System.out.print("\nThe algorithm  LSPG is starting...\nComputing...\n");
	  
	  System.out.println("Number of prototypes, result set = "+numberOfPrototypes+ "\n");
	  
	  if(numberOfPrototypes < trainingDataSet.getPosibleValuesOfOutput().size()){
		  System.out.println("Number of prototypes less than the number of clases");
		  numberOfPrototypes = trainingDataSet.getPosibleValuesOfOutput().size();
	  }
	  
	  /*
	   * Parameters:
	   * 	perc: [0,100] - percentage of reduction from original set: 
	   *    rho: [0,0.4] - constant value that indicates the movement on an attribute; it assumed the solution attributes are in [0,1].
	   *    mode: {cap, toroidal,reflection} - when creating a perturbation and goes over the range.. cap will simple not allow greater or less than 0,1; toroidal will go around the range  and relfection will change the sign of the perturbations.
	   *    MaxEval should be set to 5000*N..
	   * */
	  
	  // Step 0: random initialisation of a potential selection of prototypes from TR; 
	  // size of the solution depends on the input parameter `numberOfPrototypes`.
	  
	  PrototypeSet solution = selecRandomSet(numberOfPrototypes,true).clone() ;
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
	  
	  while (totalEval < this.MaxEval & this.rho > 0.0001){
		  
		  int numEvals = 0;
		  double fitness = fitnessDelta(solution, PublicVar.NONMODIFIED);  // no need to format as nomimal Population as they are originally selected from TR.
		  //double fitness_check2 = fitnessComputation(solution);
		  //System.err.println("\n% de acierto en training Nominal " + KNN.classficationAccuracy(solution,trainingDataSet,1)*100./trainingDataSet.size() );
			
		  if(totalEval==0)
			  System.out.println(totalEval+"\t"+fitness);

		  if(fitness == 100){
			  System.out.println("\n ****** \n optimal value for the fitness");
			  break;
		  }
		  
		  numEvals++;
		  int sizeVector = solution.size()*solution.get(0).numberOfInputs();
		 // System.out.println("Size of the vector: "+sizeVector+ "\n");
		  
		  // version 1:  iterative... from instance 0 to instance N
		  // version 2: randomly selected instances. 
		  int randomIndexes [] = UtilFunc.createShuffledIndexArr(sizeVector);
		  int goodMoves=0;
		  double currentFitness = fitness;
		  for (int iter=0;iter<sizeVector; iter++){	  
			  this.Reflection = false;

			  int pos = (suffleIndex == false) ? iter : randomIndexes[iter]; //Select the position using a random or iterative
			  //int pos =  iter;
			  
			  // calculate prototype and feature.
			  int proto = pos/numFeatures;
			  int feature= pos%numFeatures;
			  instanceIdx = pos / solution.get(0).numberOfInputs();
			  double currentValue = solution.get(proto).getInput(feature);
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
			  if (numEvals == 2) {
				  int change = 1 ;
			  }
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
				  fitCollect.add_Eval_Fit(totalEval+numEvals, fitness);
			  }else{
				  // we are going to try another correction:
				  if(!Reflection)
					  newValue = currentValue + this.rho/2.0;
				  else
					  newValue = currentValue - this.rho/2.0;
				  
				  // Check boundaries:
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
					  fitCollect.add_Eval_Fit(totalEval+numEvals, fitness);

				  }else{
					  solution.get(proto).setInput(feature,currentValue); // restore value
					  //need to update the distMatrix, reverse to the state when the feature was not changed. 
				      for (int i=0; i<trainingDataSet.size(); i++)
				      {
				    	  distMatrix[i][instanceIdx] =backup_col[i];
				      }
				      
				  }
	
			  }
			 
		  }
		  
		  if(currentFitness == fitness){
			  sameFitness++;
		  }
		  // Step 5: check if we have done something with solution.
		  if (goodMoves ==0 | sameFitness > 10){
			  sameFitness =0;
		  	  this.rho/=2;
		  }
		  totalEval +=  numEvals;
	  }// end  main-while loop
	  
	  int evalSaved = (int)(this.MaxEval-totalEval>0? this.MaxEval-totalEval: 0);
	   fitCollect.collectEvalSaved(this.MaxEval, evalSaved); //Because while loop can stop with rho conditions, 
	   //We may save a number of (this.MaxEval-totalEval) evaluations.
        System.out.println("Finish - LSPG");
        
        System.out.println("I should have done a max evaluations: "+this.MaxEval);
        System.out.println("Total number of evaluations USED: "+totalEval);
        System.out.println("Total number of evaluations SAVED: "+evalSaved);
		PrototypeSet nominalPopulation = new PrototypeSet();
        nominalPopulation.formatear(solution);
		System.err.println("\n% de acierto en training Nominal " + KNN.classficationAccuracy(nominalPopulation,trainingDataSet,1)*100./trainingDataSet.size() );
				
		return nominalPopulation;
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
      Parameters.setUse("LSPG", "<seed> <Number of neighbors>\n<Swarm size>\n<Particle Size>\n<MaxIter>\n<DistanceFunction>");        
      Parameters.assertBasicArgs(args);
      
      PrototypeSet training = PrototypeGenerationAlgorithm.readPrototypeSet(args[0]);
      PrototypeSet test = PrototypeGenerationAlgorithm.readPrototypeSet(args[1]);
      
      
      long seed = Parameters.assertExtendedArgAsInt(args,2,"seed",0,Long.MAX_VALUE);
      LSPG_DTGenerator.setSeed(seed);
      
      int k = Parameters.assertExtendedArgAsInt(args,3,"number of neighbors", 1, Integer.MAX_VALUE);
      String mode = Parameters.assertExtendedArg(args,4,"mode", 1, Integer.MAX_VALUE);
      double rho = Parameters.assertExtendedArgAsInt(args,5,"c1", 1, Double.MAX_VALUE);
      int perc = Parameters.assertExtendedArgAsInt(args,6,"particle size", 1, Integer.MAX_VALUE);
      int maxEval = Parameters.assertExtendedArgAsInt(args,7,"max iter", 1, Integer.MAX_VALUE);


      
      LSPG_DTGenerator generator = new LSPG_DTGenerator(training, k,perc, mode, rho,maxEval);
      
  	  
      PrototypeSet resultingSet = generator.execute();
      
  	//resultingSet.save(args[1]);
      //int accuracyKNN = KNN.classficationAccuracy(resultingSet, test, k);
      int accuracy1NN = KNN.classficationAccuracy(resultingSet, test);
      generator.showResultsOfAccuracy(Parameters.getFileName(), accuracy1NN, test);
  }

}
