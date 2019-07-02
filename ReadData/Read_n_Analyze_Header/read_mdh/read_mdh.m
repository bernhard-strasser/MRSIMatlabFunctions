%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%    FUNCTION TO READ IN THE header of .dat IMAGING FILES    %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% fredir = frequency encoding direction; phadir = phase encoding direction



function sMDH = read_mdh(file)




%% 0. Preparations


% READ SIZE OF HEADER
%file = fopen(sprintf('%s', image_dat_file),'r');
%headersize = fread(file,1, 'int32');


% READ CENTER OF K-SPACE, MEASURED POINTS IN FREQ ENCOD AND PHASE ENCOD DIRECTION
%fseek(file, headersize,'bof');

% DMA_length = fread(file,1,'ubit25')
% pack_bit = fread(file,1,'ubit1')
% 
% 
% chak_header = fread(file, 39, 'uint16')
% k_fredir_center = chak_header(33)+1;
% fredir_measured = chak_header(15);
% total_channel_no = chak_header(16);
% phadir_measured = chak_header(39)*2;




% phadir_measured = 1;
% fseek(file, headersize,'bof');
% while(~feof(file))
%     
%         chak_header = fread(file, 64, 'uint16')
%         if(chak_header(39) > phadir_measured)
%             phadir_measured = chak_header(39)
%         end
%         fseek(file,fredir_measured*2*4*total_channel_no + 0*(total_channel_no-1)*128,'cof');
%               
% end


% % READ TOTAL CHANNEL NUMBER
% fseek(file, -256,'eof');
% chak_header = fread(file, 16, 'uint16');
% total_channel_no = chak_header(16);






% % typedef struct
% % {
% %   unsigned short  ushLine;                  /* line index                   */
% %   unsigned short  ushAcquisition;           /* acquisition index            */
% %   unsigned short  ushSlice;                 /* slice index                  */
% %   unsigned short  ushPartition;             /* partition index              */
% %   unsigned short  ushEcho;                  /* echo index                   */	
% %   unsigned short  ushPhase;                 /* phase index                  */
% %   unsigned short  ushRepetition;            /* measurement repeat index     */
% %   unsigned short  ushSet;                   /* set index                    */
% %   unsigned short  ushSeg;                   /* segment index  (for TSE)     */
% %   unsigned short  ushIda;                   /* IceDimension a index         */
% %   unsigned short  ushIdb;                   /* IceDimension b index         */
% %   unsigned short  ushIdc;                   /* IceDimension c index         */
% %   unsigned short  ushIdd;                   /* IceDimension d index         */
% %   unsigned short  ushIde;                   /* IceDimension e index         */
% % } sLoopCounter;                             /* sizeof : 28 byte             */
% % 
% % typedef struct
% % {
% %   float  flSag;
% %   float  flCor;
% %   float  flTra;
% % } sVector;
% % 
% % 
% % typedef struct
% % {
% %   sVector  sSlicePosVec;                    /* slice position vector        */
% %   float    aflQuaternion[4];                /* rotation matrix as quaternion*/
% % } sSliceData;                               /* sizeof : 28 byte             */
% % /*--------------------------------------------------------------------------*/
% % /*  Definition of cut-off data                                              */
% % /*--------------------------------------------------------------------------*/
% % typedef struct
% % {
% %   unsigned short  ushPre;               /* write ushPre zeros at line start */
% %   unsigned short  ushPost;              /* write ushPost zeros at line end  */
% % } sCutOffData;
% % 
% % 
% % 
% % /*--------------------------------------------------------------------------*/
% % /*  Definition of measurement data header                                   */
% % /*--------------------------------------------------------------------------*/
% % typedef struct
% % {
% %   unsigned long  ulDMALength;                  // DMA length [bytes] must be                        4 byte                                               // first parameter                        
% %   long           lMeasUID;                     // measurement user ID                               4     
% %   unsigned long  ulScanCounter;                // scan counter [1...]                               4
% %   unsigned long  ulTimeStamp;                  // time stamp [2.5 ms ticks since 00:00]             4
% %   unsigned long  ulPMUTimeStamp;               // PMU time stamp [2.5 ms ticks since last trigger]  4
% %   unsigned long  aulEvalInfoMask[MDH_NUMBEROFEVALINFOMASK]; // evaluation info mask field           8
% %   unsigned short ushSamplesInScan;             // # of samples acquired in scan                     2
% %   unsigned short ushUsedChannels;              // # of channels used in scan                        2   =32
% %   sLoopCounter   sLC;                          // loop counters                                    28   =60
% %   sCutOffData    sCutOff;                      // cut-off values                                    4           
% %   unsigned short ushKSpaceCentreColumn;        // centre of echo                                    2
% %   unsigned short ushDummy;                     // for swapping                                      2
% %   float          fReadOutOffcentre;            // ReadOut offcenter value                           4
% %   unsigned long  ulTimeSinceLastRF;            // Sequence time stamp since last RF pulse           4
% %   unsigned short ushKSpaceCentreLineNo;        // number of K-space centre line                     2
% %   unsigned short ushKSpaceCentrePartitionNo;   // number of K-space centre partition                2
% %   unsigned short aushIceProgramPara[MDH_NUMBEROFICEPROGRAMPARA]; // free parameter for IceProgram   8   =88
% %   unsigned short aushFreePara[MDH_FREEHDRPARA];// free parameter                          4 * 2 =   8   
% %   sSliceData     sSD;                          // Slice Data                                       28   =124
% %   unsigned long	 ulChannelId;                  // channel Id must be the last parameter             4
% %   } sMDH; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sMDH.ulDMALength=fread(file,1,'uint32');
sMDH.lMeasUID=fread(file,1,'int32');
sMDH.ulScanCounter=fread(file,1,'uint32');
sMDH.ulTimeStamp=fread(file,1,'uint32');
sMDH.ulPMUTimeStamp=fread(file,1,'uint32');
sMDH.aulEvalInfoMask=fread(file,2,'uint32');    %8-byte, 2 elements
sMDH.ushSamplesInScan=fread(file,1,'uint16');
sMDH.ushUsedChannels=fread(file,1,'uint16');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sMDH.sLC.ushLine=fread(file,1,'uint16');                %28-byte loop counter
sMDH.sLC.ushAcquisition=fread(file,1,'uint16');
sMDH.sLC.ushSlice=fread(file,1,'uint16');
sMDH.sLC.ushPartition=fread(file,1,'uint16');
sMDH.sLC.ushEcho=fread(file,1,'uint16');
sMDH.sLC.ushPhase=fread(file,1,'uint16');
sMDH.sLC.ushRepetition=fread(file,1,'uint16');
sMDH.sLC.ushSet=fread(file,1,'uint16');
sMDH.sLC.ushSeg=fread(file,1,'uint16');
sMDH.sLC.ushIda=fread(file,1,'uint16');
sMDH.sLC.ushIdb=fread(file,1,'uint16');
sMDH.sLC.ushIdc=fread(file,1,'uint16');
sMDH.sLC.ushIdd=fread(file,1,'uint16');
sMDH.sLC.ushIde=fread(file,1,'uint16');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sMDH.sCutOff.ushPre=fread(file,1,'uint16');
sMDH.sCutOff.ushPost=fread(file,1,'uint16');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sMDH.ushKSpaceCentreColumn=fread(file,1,'uint16');
sMDH.ushDummy=fread(file,1,'uint16');
sMDH.fReadOutOffcentre=fread(file,1,'float');
sMDH.ulTimeSinceLastRF=fread(file,1,'uint32');
sMDH.ushKSpaceCentreLineNo=fread(file,1,'uint16');
sMDH.ushKSpaceCentrePartitionNo=fread(file,1,'uint16');
sMDH.aushIceProgramPara=fread(file,4,'uint16');
sMDH.aushFreePara=fread(file,4,'uint16');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sMDH.sSD.sSlicePosVec.flSag=fread(file,1,'float');
sMDH.sSD.sSlicePosVec.flCor=fread(file,1,'float');
sMDH.sSD.sSlicePosVec.flTra=fread(file,1,'float');
sMDH.sSD.aflQuaternion=fread(file,4,'float');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sMDH.ulChannelId=fread(file,1,'uint32');
sMDH.ulChannelId=fread(file,1,'uint16');
sMDH.ulDummyChannelId=fread(file,1,'uint16');
























%fclose(file);

