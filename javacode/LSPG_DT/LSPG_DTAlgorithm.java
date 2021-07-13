//
//  Main.java
//
//  Salvador García López
//
//  Created by Salvador García López 10-7-2004.
//  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
//

package keel.Algorithms.Instance_Generation.LSPG_DT;

import keel.Algorithms.Instance_Generation.Basic.PrototypeSet;
import keel.Algorithms.Instance_Generation.Basic.PrototypeGenerationAlgorithm;
import keel.Algorithms.Instance_Generation.*;
import keel.Algorithms.Instance_Generation.utilities.*;

import java.util.*;

/**
 * LSPG algorithm calling.
 * @author Isaac Triguero
 */
public class LSPG_DTAlgorithm extends PrototypeGenerationAlgorithm<LSPG_DTGenerator>
{
    LSPG_DTAlgorithm() {
		super();
		
	}

	/**
     * Builds a new ChenGenerator.
     * @param train Training data set.
     * @param params Parameters of the method.
     */
    protected LSPG_DTGenerator buildNewPrototypeGenerator(PrototypeSet train, Parameters params)
    {
       return new LSPG_DTGenerator(train, params);   
       
    }
    
     /**
     * Main method. Executes LSPG algorithm.
     * @param args Console arguments of the method.
     */
   
    public static void main(String args[]) {
    	LSPG_DTAlgorithm lspg = new LSPG_DTAlgorithm();
    	lspg.execute(args);
    	lspg.computeRuntime();
        
    }
}
