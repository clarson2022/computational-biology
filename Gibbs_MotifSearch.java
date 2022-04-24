import java.util.*;

public class Gibbs_MotifSearch extends EM_MotifSearch{
    
    /**
     * constructor
     */
    public Gibbs_MotifSearch(java.lang.String fileName, int motifLength){
        super(fileName, motifLength);
    }

    /**
     * The Expectation step in the EM algorithm.
     * Based on the matrix model, identifies motif instances in the sequences. One motif instance is 
     * identified in each sequence. For each sequence, a motif instance is chosen by sampling the scores 
     * of each possible motif instance in that sequence. The score of each possible motif instance is 
     * based on the matrix model.
     */
    public void determineMotifInstances(){
        //for each sequence
        for(int i = 0; i<numSequences; i++){
            //choose motif sequence
            Vector<Double> scores = getScores(sequences.get(i));
            //use this vector to pass to getIndexViaSampling
            int index = getIndexViaSampling(scores);
            //the index this returns is what is passed to 
            instanceLocations.set(i, index);
        }
    }
       
    private Vector<Double> getScores(String s){
        //get the possible scores of a sequence and put them in vector
        Vector<Double> scores = new Vector<Double>();
        for(int i = 0; i<s.length()-motifLength+1; i++){
            String current = s.substring(i, i+motifLength);
            scores.add(getScoreForMotifInstance(current));
        }
        return scores;
    }
    

    /**
     * Returns the index of a randomly sampled value in a Vector.
     * One value from the Vector is chosen at random and the value's index (not the value 
     * itself) is returned.
     */
    public int getIndexViaSampling(Vector<Double> values){
        //normalize values by summing and dividing each entry by the sum
        Vector<Double> normalized = normalizeValues(values);
        //Convert the values from a probability distribution to a cumulative distribution.
        Vector<Double> cumulative = convertToCumulative(normalized);
        //Generate a number uniformly at random between 0.0 and 1.0. Return the index of the smallest value in the cumulative distribution that is at least as big as the random number.
        int index = selectIndex(cumulative);
        return index;
    }

    private static Vector<Double> normalizeValues(Vector<Double> values){
        double sum = 0;
        double normalizedVal; 
        for(double num:values){
            sum+= num;
        }
        for(int i = 0; i<values.size(); i++){
            normalizedVal = values.get(i)/sum;
            values.set(i, normalizedVal);
        }
        return values;
    }

    private static Vector<Double> convertToCumulative(Vector<Double> values){
        Vector<Double> cumulative = new Vector<Double>();
        double sumSoFar = 0;
        //go through each in value and add the sum to far
        for(double num:values){
            sumSoFar+=num;
            cumulative.add(sumSoFar);
        }
        return cumulative;
    }

    private static int selectIndex(Vector<Double> cumulativeDist){
        //generate number btwn 0 and 1
        Random rand = new Random();
        double random = rand.nextDouble();
        int index = 0;
        for(int i = 0; i<cumulativeDist.size(); i++){
            if(cumulativeDist.get(i)>=random){
                index = i;
                break;
            }
        }
        return index;
    }


    public static void main(String[] args) {
        Gibbs_MotifSearch g = new Gibbs_MotifSearch(args[0], Integer.parseInt(args[1]));
        
        g.run_EM_multiple_times(Integer.parseInt(args[2]));

        //The set of motif instances, matrix, consensus sequence, and information content corresponding to the maximum over all iterations are output.
        
        System.out.println(g.motifInstancesToString());
        System.out.println(g.matrixToString());
        System.out.println(g.getConsensusSequence());
        System.out.println(g.getInformationContentOfMatrix());   
    }
}