
/**
 * Hierarchical_Clustering executes the hierarchical clustering algorithm. 
 *
 * @author Catherine Larson and Skylar Kolisko
 * @version April 7 2022
 */
public class Hierarchical_Clustering extends Clustering
{
    // instance variables - replace the example below with your own
    int numClusters; 
    
    /**
     * Constructor for objects of class Hierarchical_Clustering
     */
    public Hierarchical_Clustering(String fileName, int numClust)
    {
        // initialise instance variables
        super(fileName);
        numClusters = numClust;
        hierarchical();
    }

    public void hierarchical(){
        initiallyAssignOneGeneToEachCluster();
        while(numClusters!=clusters.size()){
            mergeTwoClosestClusters();
        }
    }
    
    public void initiallyAssignOneGeneToEachCluster(){
        for(int i = 0; i<genes.size(); i++){
            //make new cluster to add to clusters vector
            Cluster newCluster = new Cluster();
            //add gene to it's own cluster 
            newCluster.addGene(genes.get(i));
            //add cluster to vector of all clusters
            clusters.add(newCluster);
            System.out.println(newCluster);
        }
    }
    
    private double findDistance(Cluster cluster1, Cluster cluster2){
        return cluster1.getDistanceToCluster(cluster2);
    }
    
    public void mergeTwoClosestClusters(){
        double smallest = -1;
        double distance;
        int cA = -1; 
        int cB = -1; 
        for(int i = 0; i<clusters.size(); i++){
            if (smallest ==-1){
                distance = findDistance(clusters.get(i), clusters.get(i+1)); 
                smallest = distance;
                cA = i;
                cB = i+1;
            }
            for (int j = i+1; j<clusters.size(); j++){
                distance = findDistance(clusters.get(i), clusters.get(j));
                if (distance<smallest){
                    smallest = distance;
                    cA = i;
                    cB = j;
                }
            }
        }        
        clusters.get(cA).absorbCluster(clusters.get(cB));
        clusters.remove(cB);
    }
    
    public static void main(String[] args){
        int numClusters = Integer.parseInt(args[1]);
        Hierarchical_Clustering h = new Hierarchical_Clustering(args[0], numClusters);
        System.out.println(h.toString());
    }
    
}
