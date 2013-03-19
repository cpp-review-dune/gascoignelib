/*----------------------------   tsuchimikado.h     ---------------------------*/
/*      $Id: tsuchimikado.h,v 1.2 2008/10/13 11:44:28 richter Exp $                 */
#ifndef __tsuchimikado_H
#define __tsuchimikado_H
/*----------------------------   tsuchimikado.h     ---------------------------*/



namespace Tsuchimikado
{

  /**
   *
   * see description later on for explanation
   *
   **/
  
  const int STD_BOUNDARY_QUAD_DIRECTION[6][4]={{0,1,2,3},
					       {1,5,6,2},
					       {2,6,7,3},
					       {0,4,5,1},
					       {4,7,6,5},
					       {0,3,7,4}};

  static const int MESH_IS_ISOTROPIC = 1;
}



	  /**
	   *
	   *    LINE
	   *
	   *  0-------1
	   *
	   *    REFINEMENT
	   *
	   * REF_X
	   *
	   *    0   1
	   *  0---x---1
	   *
	   *  type=0 no children
	   *  type=1 refined, 2 children
	   *
	   **/

	  /**
	   *
	   *    QUAD
	   *
	   * 3-------2
	   * |   2   |
	   * |3     1|
	   * |   0   |
	   * 0-------1
	   *
	   *    REFINEMENT
	   *
	   * REF_X:     REF_Y:     REF_X&REF_Y   
	   *		                    
	   * 3-------2	3-------2  3-------2
	   * |   |   |	|   1   |  | 3 | 2 |
	   * | 0 | 1 |	|-------|  |---|---|
	   * |   |   |	|   0   |  | 0 | 1 |
	   * 0-------1	0-------1  0-------1
	   *
	   *  type=1      type=2     type=3
	   *
	   * refined edges have the same orientation as the father,
	   *
	   *  for all types: children sorted according to node with
	   *  lowest number
	   *
	   **/

	  /**
	   *    HEX
	   *       7-------6        x-------x        x-------x 
	   *      /.      /|       /       /|       /|       | 
	   *     / .     / |      /       / |      / |   4   | 
	   *    /  .    /  |     /   2   /  |     /  |       | 
	   *   /   4.../...5    /       /   x    /   x-------x 
	   *  /   .   /   /    /       / 1 /    / 5 /       /  
	   * 3-------2   /    x-------x   /    x   /       /   
	   * |  .    |  /     |       |  /     |  /   3   /    
	   * | .     | /      |   0   | /      | /       /     
	   * |.      |/	      |       |/       |/       /      
	   * 0-------1        x-------x        x-------x       
	   *
	   * QUAD standard ordering (not a must):
	   * viewed from the outer side: against clockwise,
	   * lowest number first.
	   *
	   *  0    0 1 2 3
	   *  1    1 5 6 2
	   *  2    2 6 7 3
	   *  3    0 4 5 1 
	   *  4    4 7 6 5
	   *  5    0 3 7 4
	   *
	   * LINE ordering
	   *
	   *       .---6---.  
	   *      /.      /|  
	   *     / 7     / 5  
	   *    10  .   9  |  
	   *   /   ...4/....  
	   *  /   .   /   /   
	   * .---2---.   /    
	   * | 11    |  8     
	   * 3 .     1 /      
	   * |.      |/	      
	   * .---0---.        
	   *
	   *
	   *
	   *    REFINEMENT
	   *
	   *       REF_X:           REF_Y:           REF_Z:    
	   *		                                       
	   *       x-------x        x-------x        x-------x 
	   *      /   /   /|       /       /|       /       /| 
	   *     /   /   / |      /       / |      /   1   / | 
	   *    /   /   /  |     /       / /|     /-------/  | 
	   *   /   /   /   x    /       / / x    /       /|  x 
	   *  /   /   /   /    /       / / /    /   0   / | /  
	   * x-------x   /    x-------x / /    x-------x  |/   
	   * |   |   |  /     |   1   |/ /     |       |  /    
	   * | 0 | 1 | /      |-------| /      |       | /     
	   * |   |   |/	      |   0   |/       |       |/      
	   * x-------x        x-------x        x-------x       
	   *
	   *    REF_X&REF_Y:     REF_X&REF_Z:       REF_Y&REF_Z:    
	   *		                                       
	   *       x-------x        x-------x        x-------x 
	   *      /   /   /|       /   /   /|       /       /| 2
	   *     /   /   / |      / 3 / 2 / |      /       / | 
	   *    /   /   / /|     /-------/  |     /-------/ /| 3
	   *   /   /   / / x    /   /   /|  x    /       /|/ x 
	   *  /   /   / / /    / 0 / 1 / | /    /       / | /  
	   * x-------x / /    x-------x  |/    x-------x /|/   
	   * | 3 | 2 |/ /     |   |   |  /     |   1   |/ /    
	   * |---+---| /      |   |   | /      |-------| /     
	   * | 0 | 1 |/	      |   |   |/       |   0   |/      
	   * x-------x        x-------x        x-------x       
	   *
	   *
	   *  REF_X&REF_Y&REF_Z: same numbering as nodes
	   *
	   *  for all types: children sorted according to node with
	   *  lowest number
	   *
	   **/


	  /**
	   *
	   * MASTER / SLAVE
	   *
	   * If an Element E is refined anisotropic, the master/slave
	   * is not unique any more. It always indicates to the finer
	   * element.
	   *
	   *
	   *  ------- 	
	   * |   |   |	Element K, children K0, K1,
	   * | 0 | 1 L	master/slave of L is K1, not K
	   * |   |   |	
	   *  ------- 	
	   *
	   **/


           /**
	    *
	    * HANGING
	    *
	    *     2D:
	    * line hangs, if it is boundary of active quad and refined
	    *
	    *     3D:
	    * line hangs, if it is boundary of active hex and refined
	    *	    
	    * quad hangs, if it is boundary of active hex and refined
	    *  isotropically, or refined anisotropically and both
	    *  children are refined in the other anisotropic direction:
	    *   
	    * *-------*	
	    * |   |   |	
	    * |-- x --|	 
	    * |   |   |	
	    * *-------*
	    *
	    **/



/*----------------------------   tsuchimikado.h     ---------------------------*/
/* end of #ifndef __tsuchimikado_H */
#endif
/*----------------------------   tsuchimikado.h     ---------------------------*/
