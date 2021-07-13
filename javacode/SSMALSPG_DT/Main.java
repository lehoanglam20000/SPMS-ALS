//
//  Main.java
//
//  Salvador Garc�a L�pez
//
//  Created by Salvador Garc�a L�pez 3-10-2005.
//  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
//

package keel.Algorithms.Instance_Generation.SSMALSPG_DT;

public class Main {

  public static void main (String args[]) {

    SSMALSPG_DT ssma;

    if (args.length != 1)
      System.err.println("Error. A parameter is only needed.");
    else {
      ssma = new SSMALSPG_DT (args[0]);
      ssma.Script = args[0];
      ssma.ejecutar();
    }
  }
}
