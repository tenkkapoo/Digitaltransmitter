module tb_ht_chip #(
    g_Rs = 100000000.0,
    g_file_A = "/home/pro/hattrick/spelman/hattrick_thesydekick/Entities/ht_chip/Simulations/rtlsim/A_tmpoz7bpsy2.txt",
    g_file_Z = "/home/pro/hattrick/spelman/hattrick_thesydekick/Entities/ht_chip/Simulations/rtlsim/Z_tmp3xa5_yrq.txt",
    g_file_control_write = "/home/pro/hattrick/spelman/hattrick_thesydekick/Entities/ht_chip/Simulations/rtlsim/control_write_tmp36tb352z.txt"
);
//timescale 1ps this should probably be a global model parameter
//Parameter definitions
parameter integer c_Ts=1/(g_Rs*1e-12);
//Register definitions
reg reset;
reg A;
reg clock;
reg initdone;

//Wire definitions
wire Z;

//Assignments
//Variables for the io_files
integer status_A, f_A;
initial f_A = $fopen(g_file_A,"r");
integer status_Z, f_Z;
initial f_Z = $fopen(g_file_Z,"w");
integer status_control_write, f_control_write, ctstamp_control_write, ptstamp_control_write, tdiff_control_write;
initial ctstamp_control_write=0;
initial ptstamp_control_write=0;
integer buffer_reset;
integer buffer_initdone;
initial f_control_write = $fopen(g_file_control_write,"r");



//DUT definition
ht_chip ht_chip (
    .reset(reset),
    .A(A),
    .Z(Z)
);
//Master clock is omnipresent
always #(c_Ts/2.0) clock = !clock;

//io_out
         always @(posedge clock)
begin
if ( ~$isunknown(Z) 
&& initdone ) begin
$fwrite(f_Z, "%d\n",
    Z
);
        end
    end


//Execution with parallel fork-join and sequential begin-end sections
initial #0 begin
fork
A = 'b0;    
clock = 'b0;    


    // Sequences enabled by initdone
    $display("Ready to read"); 
    while (!$feof(f_A)) begin
   @(posedge clock)
        if ( initdone ) begin
        status_A = $fscanf(f_A, "%d\n",
    A
);
        end
    end
begin
while(!$feof(f_control_write)) begin
    tdiff_control_write = ctstamp_control_write-ptstamp_control_write;
    #tdiff_control_write begin
        ptstamp_control_write = ctstamp_control_write;
        reset = buffer_reset;
        initdone = buffer_initdone;
        status_control_write = $fscanf(f_control_write, "%d\t%d\t%d\n",
            ctstamp_control_write,
            buffer_reset,
            buffer_initdone
        );
    end
end
tdiff_control_write = ctstamp_control_write-ptstamp_control_write;
#tdiff_control_write begin
    ptstamp_control_write = ctstamp_control_write;
    reset = buffer_reset;
    initdone = buffer_initdone;
end
end

join

//Close the io_files
$fclose(f_A);
$fclose(f_Z);
$fclose(f_control_write);


$finish;
end

endmodule