//
//  Main.java
//
//  Salvador García López
//
//  Created by Salvador García López 3-10-2005.
//  Copyright (c) 2004 __MyCompanyName__. All rights reserved.
//

package keel.Algorithms.Instance_Generation.SSMALSHADE;

public class Main {

  public static void main (String args[]) {

    SSMALSHADE ssma;

    if (args.length != 1)
      System.err.println("Error. A parameter is only needed.");
    else {
      ssma = new SSMALSHADE (args[0]);
      ssma.Script = args[0];
      ssma.ejecutar();
    }
  }
}
