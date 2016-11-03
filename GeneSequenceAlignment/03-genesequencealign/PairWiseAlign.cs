using System;
using System.Collections.Generic;
using System.Text;

namespace GeneticsLab
{
    class PairWiseAlign
    {
        int MaxCharactersToAlign;
        int INDEL = 5;
        int SUB = 1;
        int MATCH = -3;

        public PairWiseAlign()
        {
            // Default is to align only 5000 characters in each sequence.
            this.MaxCharactersToAlign = 5000;
        }

        public PairWiseAlign(int len)
        {
            // Alternatively, we can use an different length; typically used with the banded option checked.
            this.MaxCharactersToAlign = len;
        }

        /// <summary>
        /// this is the function you implement.
        /// </summary>
        /// <param name="sequenceA">the first sequence</param>
        /// <param name="sequenceB">the second sequence, may have length not equal to the length of the first seq.</param>
        /// <param name="banded">true if alignment should be band limited.</param>
        /// <returns>the alignment score and the alignment (in a Result object) for sequenceA and sequenceB.  The calling function places the result in the dispay appropriately.
        /// 
        public ResultTable.Result Align_And_Extract(GeneSequence sequenceA, GeneSequence sequenceB, bool banded)
        {
            ResultTable.Result result = new ResultTable.Result();
            int score;                                                       // place your computed alignment score here
            string[] alignment = new string[2];                              // place your two computed alignments here


           
            // ********* these are placeholder assignments that you'll replace with your code  *******
            score = 0;                                                
            alignment[0] = "";
            alignment[1] = "";
            // ***************************************************************************************
            //System.Console.WriteLine(sequenceA);
            Final_Data f = null;

            if (banded)
            {
                f = banded_align(sequenceA.Sequence, sequenceB.Sequence);
            }
            else
            {
                f = normal_align(sequenceA.Sequence, sequenceB.Sequence);
            }
            

            score = f.Score;
            alignment[0] = f.Vals.A;
            alignment[1] = f.Vals.B;
            result.Update(score,alignment[0],alignment[1]);                  // bundling your results into the right object type 
            return(result);
        }

        public Final_Data normal_align(String sequenceA, String sequenceB)
        {
            /*
            * The following is the normal way of solving
            */
            int lengthA = Math.Min(sequenceA.Length, this.MaxCharactersToAlign) + 1;
            int lengthB = Math.Min(sequenceB.Length, this.MaxCharactersToAlign) + 1;
            int[,] results = new int[lengthA, lengthB];
            char[,] path = new char[lengthA, lengthB];


            for (int i = 0; i < lengthA; i++)
            {
                results[i, 0] = 0;
            }

            for (int i = 0; i < lengthB; i++)
            {
                results[0, i] = 0;
            }

            //These nested loops give the O(mn) time, I have to loop through the inner string length once for each letter in the outer string
            for (int i = 1; i < lengthA; i++)
            {
                for (int j = 1; j < lengthB; j++)
                {
                    //System.Console.WriteLine( string.Format("IN MAIN LOOP: i {0}, j {1}", i, j ) );
                    min_return vals = calculate(results[i - 1, j] + INDEL,
                        results[i, j - 1] + INDEL,
                        results[i - 1, j - 1] + diff(sequenceA[i - 1], sequenceB[j - 1])); 

                    //System.Console.Write(string.Format(",{0},", vals.Val));
                    results[i, j] = vals.Val;
                    path[i, j] = vals.Prev;
                }
            }

            StringVals ret_strings = reverseStrings(sequenceA.Substring(0, lengthA -1), sequenceB.Substring(0, lengthB-1), path);
            //StringVals ret_strings = new StringVals("1", "2");
            //System.Console.WriteLine(string.Format("A: {0} B: {1}", ret_strings.A, ret_strings.B));


            return new Final_Data(ret_strings, results[lengthA - 1, lengthB - 1]);
        }
        public Final_Data banded_align(String sequenceA, String sequenceB)
        {
            /*
            * The following is the banded way of solving
            */
            int lengthA = Math.Min(sequenceA.Length, this.MaxCharactersToAlign) + 1;
            int lengthB = Math.Min(sequenceB.Length, this.MaxCharactersToAlign) + 1;
            int[,] results = new int[lengthA, lengthB];
            char[,] path = new char[lengthA, lengthB];


            //For banded I'll only need the first 4 values as 0
            for (int i = 0; i < 4; i++)
            {
                results[i, 0] = 0;
            }

            for (int i = 0; i < 4; i++)
            {
                results[0, i] = 0;
            }

            //This will be used to calculate how far to go on each row, just needs to increase by 1 each time and account for the right edge at the end
            int base_num = 5;
            //This number will be used to indicate the starting index, should just go up by 1 each time
            int start = 1;
            //This is where we get the O(n+m) time, it's because I'm not looping through the entire inner string for each value in the outer one, I'm only looping through a small amount
            for (int i = 1; i < lengthA; i++)
            {
                for (int j = start; j < base_num; j++)
                {
                    int? l = 0;
                    int? u = 0;
                    //I need to account for the unfilled parts in the array when I'm grabbing min values, only matters for left and up at the edges
                    if (j == base_num - 1)
                    {
                        //this will be a case where the up can't be grabbed
                        //I chose 10 arbitrarily so it will never be the min
                        u = null;

                    }
                    else
                    {
                        //normal behavior
                        u = results[i, j - 1] + INDEL;
                    }
                    if (j == base_num - 7)
                    {
                        //I believe this will be the case where the left can't be grabbed
                        //I chose 10 arbitrarily so it will never be the min
                        l = null;
                    }
                    else
                    {
                        //normal behavior
                        l = results[i - 1, j] + INDEL;
                    }
                    min_return vals = calculate(l,
                        u,
                        results[i - 1, j - 1] + diff(sequenceA[i - 1], sequenceB[j - 1])); //use -1 for the diff because of the empty string in out loops

                    //System.Console.Write(string.Format(",{0},", vals.Val));
                    results[i, j] = vals.Val;
                    path[i, j] = vals.Prev;
                }
                //this accounts for the beginning part
                //I also need to consider the end for this one, I need to make sure this will do it correctly
                if (base_num >= 8)
                {
                    start++;
                }
                //the base num increases by 1 until it equals the length of the B string then stops increasing to avoid any out of bounds errors
                if (base_num+1 >= lengthB )
                {
                    base_num = lengthB;
                }
                else
                {
                    base_num++;
                }
            }


            //StringVals ret_strings = reverseStrings(sequenceA.Sequence.Substring(0, lengthA -1), sequenceB.Sequence.Substring(0, lengthB-1), path);
            StringVals ret_strings = new StringVals("1", "2");

            return new Final_Data(ret_strings, results[lengthA - 1, lengthB - 1]);
        }

        public StringVals reverseStrings(String a, String b, char[,] prev)
        {
            // Console.WriteLine(prev[5, 5].ToString());

            //Console.WriteLine("A string substring: " + a);
            //Console.WriteLine("B string substring: " + b);
            //Console.WriteLine("The previous: ");
            bool done = false;

            String revA = "";
            String revB = "";

            int i = a.Length; //this will have to be the actual length
            int j = b.Length; //this will have to be the actual length

            while (false == done)
            {
                char cur = prev[i,j];

                if ( cur == 'd')
                {
                    revA += a[i - 1];
                    revB += b[j - 1];

                    i--;
                    j--;
                }
                else if (cur == 'l')
                {
                    revA += a[i - 1];
                    revB += '-';

                    i--;
                }
                else //cur is u
                {
                    revA += '-';
                    revB += b[j - 1];

                    j--;
                }

                //I'll need to check if these values are right
                //look at size of prev array vs the size of the strings
                if (i == 0 && j ==0)
                {
                    done = true;
                }
            }


            //  System.Console.WriteLine(revA.Length);
            //System.Console.WriteLine(revA);
            // System.Console.WriteLine(revB);
            /*
            string output = "";
            for (int k = revA.Length - 1; k > -1; k--)
            {
                output += revA[k];
            }

            string output2 = "";
            for (int k = revB.Length - 1; k > -1; k--)
            {
                output2 += revB[k];
            }
           
            // System.Console.WriteLine(output);
            //System.Console.WriteLine(output2);

            StringVals c = new StringVals(revA, revB);

            return c;
             */


            //Need to reverse the strings now

            /*
             char[] ra = revA.ToCharArray();
             char[] rb = revB.ToCharArray();

            Array.Reverse(ra);
             Array.Reverse(rb);
            return new StringVals(new string(ra), new string(rb));
            */
            return new StringVals(revA, revB);



        }

        public min_return calculate(int? left, int? up, int diagonal)
        {
            //if (debug) System.Console.WriteLine(string.Format("IN CALCULATE: left {0}, right {1}, diagonal {2}", left,right,diagonal  ));
            if (left != null && up != null)
            {
                if (left > up)
                {
                    if (up > diagonal)
                    {
                        //diagonal is minimun
                        //add previous value
                        return new min_return('d', diagonal);
                    }
                    else
                    {
                        return new min_return('u', (int)up);
                    }

                }
                if (up > left)
                {
                    if (left > diagonal)
                    {
                        //diagonal is minimun
                        //add previous value
                        return new min_return('d', diagonal);
                    }
                    else
                    {
                        return new min_return('l', (int)left);
                    }

                }


                //left and up are equal 
                if (left > diagonal)
                {
                    //diagonal is minimun
                    //add previous value
                    return new min_return('d', diagonal);
                }

                return new min_return('l', (int)left);
            }
            if (left == null)
            {
                //only need to compare up and diagonal
                if (up > diagonal)
                {
                    //diagonal is minimun
                    //add previous value
                    return new min_return('d', diagonal);
                }
                else
                {
                    return new min_return('u', (int)up);
                }

            }
            if (up == null)
            {
                //only need to compare left and diagonal
                if (left > diagonal)
                {
                    //diagonal is minimun
                    //add previous value
                    return new min_return('d', diagonal);
                }
                else
                {
                    return new min_return('l', (int)left);
                }
            }
            return null;
            
        }

        /*
         *  Computes the diff and returns the specified values based on result
         */
        public int diff( Char a, Char b )
        {
            if (0 == a.CompareTo(b))
            {
                return MATCH;
            }
            else
            {
                return SUB;
            }
            
        }
    }

    /*
     * Class to hold backpointer value and the cost
     */
    public class min_return
    {
        char prev;
        int val;
        public min_return(char prev, int val)
        {
            this.Prev = prev;
            this.val = val;
        }

        public char Prev
        {
            get
            {
                return prev;
            }

            set
            {
                prev = value;
            }
        }

        public int Val
        {
            get
            {
                return val;
            }

            set
            {
                val = value;
            }
        }
    }

    /*
     * Class to hold two strings
     */
    public class StringVals
    {
        String a;
        String b;
        public StringVals(String a, String b)
        {
            this.A = a;
            this.B = b;
        }

        public string A
        {
            get
            {
                return a;
            }

            set
            {
                a = value;
            }
        }

        public string B
        {
            get
            {
                return b;
            }

            set
            {
                b = value;
            }
        }
    }


    /*
     * Class to hold StringVals and Score
     */
    public class Final_Data
    {
        StringVals vals;
        int score;
        public Final_Data(StringVals a, int s)
        {
            this.Vals = a;
            this.Score = s;
        }

        public StringVals Vals
        {
            get
            {
                return vals;
            }

            set
            {
                vals = value;
            }
        }

        public int Score
        {
            get
            {
                return score;
            }

            set
            {
                score = value;
            }
        }

        
        
    }
}
