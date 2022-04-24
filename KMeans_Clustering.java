import java.util.*;
/**
 * KMeans_Clustering executes the k means clustering algorithm. 
 *
 * @author Catherine Larson and Skylar Kolisko
 * @version April 7 2022
 */
public class KMeans_Clustering extends Clustering
{
    int k;
    Vector<Vector<Double>> means;
    double oldDistances;
    /**
     * Constructor for objects of class KMeans_Clustering
     */
    public KMeans_Clustering(String fileName, int numClust)
    {
        // initialise instance variables
        super(fileName);
        k = numClust;
        means = new  Vector<Vector<Double>>();
        kMeans();
    }

    public void initializeAllClusters(){
        for(int i=0; i<clusters.size(); i++){
            clusters.get(i).initialize();
        }
    }
    
    public void randomlyAssignGenesToClusters(){
        Vector<Integer> choices = new Vector<Integer>();
        for(int i = 0; i<genes.size(); i++){
            choices.add(i);
        }
        Random rand = new Random();  
        Gene currentGene;
        int counter = 0;
        while(choices.size()>0){
            int currentCluster = counter % k;
            
            int randomInt = rand.nextInt(choices.size());
            
            int geneToGet = choices.get(randomInt);
            currentGene = genes.get(geneToGet);
            clusters.get(currentCluster).addGene(currentGene);
                    
            choices.remove(randomInt);
            counter++;
        }
        
    }
    
    public Vector<Vector<Double>> getMeansOfAllClusters(){
        Vector<Vector<Double>> means = new  Vector<Vector<Double>> ();
        //calculate
        //for each cluster
        for(int i= 0; i< clusters.size(); i++){
            //getClusterMean
            //add it to means
            means.add(clusters.get(i).getClusterMean());
        }
        return means;
    }
    
    public boolean assignGenesToClusters(){        
        //clear out all clusters
        initializeAllClusters();
        double dist;
        int clusterIndex = 0;
        double totalDist = 0;
        boolean done_improving=false;
        //for each gene
        for(int i = 0; i < genes.size(); i++){
            //assign to a cluster
            double lowest = genes.get(i).distanceToExpressionVector(means.get(0));
            for(int j = 1; j< means.size(); j++){
                dist = genes.get(i).distanceToExpressionVector(means.get(j));
                if(dist<lowest){
                    lowest = dist;
                    clusterIndex = j;
                }
            }
            clusters.get(clusterIndex).addGene(genes.get(i));
            totalDist+=lowest;
        }
        
        //compare old clusters to new clusters, if better then return right boolean
        if (totalDist==oldDistances){
            done_improving=true;
        }
        else{
            oldDistances = totalDist;
        }
        
        return done_improving;
    }
        
    public void populateEmptyClusters(){
        for(int i = 0; i<clusters.size(); i++){
            if(clusters.get(i).getSizeOfCluster() == 0){
                //put a gene in it
                populate(i);
            }
        }
    }
    
    private void populate(int emptyClusterIndex){
        //in case of invalid input
        Random rand = new Random();
        //int count = 0;
        int randomCluster = rand.nextInt(clusters.size());
        while(clusters.get(randomCluster).getSizeOfCluster()<2){
            randomCluster = rand.nextInt(clusters.size());
        }
        Cluster randClust = clusters.get(randomCluster);
        int randomGeneIndex = rand.nextInt(randClust.getSizeOfCluster());
        Gene  replacementGene = randClust.getGene(randomGeneIndex);
        
        clusters.get(emptyClusterIndex).addGene(replacementGene);
        
        randClust.removeGene(replacementGene);
        
    }
    
    public void kMeans(){
        for(int i=0; i<k; i++){
            Cluster clusterK = new Cluster();
            clusters.add(clusterK);
        }        
        randomlyAssignGenesToClusters();
        boolean converged = false;
        
        while(!converged){
            //calculate means of each cluster
            means = getMeansOfAllClusters();            
            // assign each gene to cluster with closest mean
            converged = assignGenesToClusters();
            //iff assign returns true then converged = true
            // check for empty clusters
            populateEmptyClusters();
        }
    }
    
    public static void main(String[] args){
        int k = Integer.parseInt(args[1]);
        KMeans_Clustering c = new KMeans_Clustering(args[0], k);
        System.out.println(c.toString());
    }
    
}
