import java.util.*;
/**
 * Hierarchical_Clustering executes the hierarchical clustering algorithm. 
 *
 * @author Catherine Larson and Skylar Kolisko
 * @version April 7 2022
 */
public class CAST_Clustering extends Clustering
{
    double t;
    HashMap<Gene, Integer> unassignedGenes = new HashMap<Gene, Integer>();
    
    /**
     * Constructor for objects of class CAST_Clustering
     */
    public CAST_Clustering(String fileName, double threshold)
    {
        super(fileName);
        t = threshold;
        //all genes start unassigned
        for(int i = 0; i<genes.size(); i++){
            unassignedGenes.put(genes.get(i), 1);
        }
        cast();
    }
    
    public int addGenesWithHighAffinity(Cluster current){
        int genesAdded = 0;
        //get the cluster average
        Vector<Double> clusterMean = current.getClusterMean();
        //go through each gene  
        for(int i = 0; i < genes.size(); i++){
            Gene currentGene = genes.get(i);
            //check if it's in any cluster
            if(unassignedGenes.containsKey(currentGene)){
                //if it's close to the mean
                if(currentGene.distanceToExpressionVector(clusterMean) <= t){
                    //then add the gene and increment count
                    current.addGene(currentGene);
                    unassignedGenes.remove(currentGene);
                    genesAdded++;
                }
            }
        }
        return genesAdded;
    }
    
    public int removeGenesWithLowAffinity(Cluster current){
        Vector<Double> clusterMean = current.getClusterMean(); //recalculate mean here, once
        int genesRemoved = 0;
        //for every gene in the cluster
        for(int i = 0; i<current.getSizeOfCluster(); i++){
            //check if it's still within the range
            Gene currentGene = current.getGene(i);
            if(currentGene.distanceToExpressionVector(clusterMean)>t){
                current.removeGene(currentGene);
                unassignedGenes.put(currentGene, 1);
                genesRemoved ++;
            }
        }  
        return genesRemoved;
    }
    
    private Gene getUnassignedGene(){
        Gene selected = genes.get(0); 
        boolean unselected = true;
        int i = 0;
        while(unselected){
            selected = genes.get(i);
            if(unassignedGenes.containsKey(selected)){
                unselected = false;
            }
            i++;
        }
        return selected;
    }
    
    public void cast(){
        int genesAdded, genesRemoved;
        boolean clusterUnfinished = true;
        
        while(unassignedGenes.size()>0){
            //get a new gene that's not in a cluster
            Gene nextGene = getUnassignedGene();
            //make a new cluster from it
            Cluster currentCluster = new Cluster();
            currentCluster.addGene(nextGene);
            clusters.add(currentCluster);
            //and finish it's cluster
            clusterUnfinished = true;
            while(clusterUnfinished){
                //run highAffinity to get everything close enough to it
                genesAdded = addGenesWithHighAffinity(currentCluster);
                //run lowAffinity
                genesRemoved = removeGenesWithLowAffinity(currentCluster);
                //if an any point both of those are not adding new things to the cluster or removing then unfinished = true
                if(genesAdded == 0 && genesRemoved == 0){
                    clusterUnfinished = false;
                }
            }
        }
    }
    
    public static void main(String[] args){
        double t = Double.parseDouble(args[1]);
        CAST_Clustering c = new CAST_Clustering(args[0], t);
        System.out.println(c.toString());
    }
}
