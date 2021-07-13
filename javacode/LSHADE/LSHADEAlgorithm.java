//
//  Main.java
//
//  Salvador García López
//
//  Created by Salvador García López 10-7-2004.
//  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
//

package keel.Algorithms.Instance_Generation.LSHADE;

import keel.Algorithms.Instance_Generation.Basic.PrototypeSet;
import keel.Algorithms.Instance_Generation.Basic.PrototypeGenerationAlgorithm;
import keel.Algorithms.Instance_Generation.*;
import keel.Algorithms.Instance_Generation.utilities.*;

import java.util.*;

/**
 * JADE algorithm calling.
 * @author Isaac Triguero
 */
public class LSHADEAlgorithm extends PrototypeGenerationAlgorithm<LSHADEGenerator>
{
    /**
     * Builds a new ChenGenerator.
     * @param train Training data set.
     * @param params Parameters of the method.
     */
    protected LSHADEGenerator buildNewPrototypeGenerator(PrototypeSet train, Parameters params)
    {
       return new LSHADEGenerator(train, params);    
    }
    
     /**
     * Main method. Executes JADE algorithm.
     * @param args Console arguments of the method.
     */
    public static void main(String args[])
    {
        LSHADEAlgorithm isaak = new LSHADEAlgorithm();
        isaak.execute(args);
    }
}
