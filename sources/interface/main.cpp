// --------------------------------------------------------------------------
// Source code provided FOR REVIEW ONLY, as part of the submission entitled
// "Moving Level-of-Detail Surfaces".
//
// A proper version of this code will be released if the paper is accepted
// with the proper licence, documentation and bug fix.
// Currently, this material has to be considered confidential and shall not
// be used outside the review process.
//
// All right reserved. The Authors
// --------------------------------------------------------------------------

#define TINYPLY_IMPLEMENTATION
#include "mainwindow.h"
#include <QApplication>

void exitUsage(const string& prog)
{
	std::cout << "--------------Help-----------------\n"
		<< "Usage: " << prog << " -i NAME_OF_INPUT_FILE [OPTIONS]\n"
		<< "Or for manual use: " << prog << "\n"
		<< "Points will be automatically projected on the input pointset file\n"
		<< "The output folder determines where the output pointset and the timings will be saved\n"
		<< "Formats supported :\n"
		<< "- ply\n"
		<< "- pn\n"
		<< "\n-------Options---------\n"
		<< "-h: display this help\n"
		<< "-o: path of output folder\n"
		<< "-f: path of output file\n"
		<< "-c: use CPU only\n"
		<< "-d: only compute the dual-contouring\n"
		<< "-a knn: use APSS with knn nearest neighbors\n"
		<< "-e r: use APSS with ball of radius r, r<0 => automatic radius\n"
		<< "-k cutoff sigma1 factor number_of_sigmas: use the multiple gaussian kernel\n"
		<< "-r cutoff exponent epsilon: use the rational kernel\n"
		<< "-s cutoff sigma: use the single gaussian kernel\n"
		<< "-g: use global APSS instead of MLoDS\n"
		<< "-m depth: adaptive octree depth\n"
		<< "-b: use octree built on input points\n"
		<< "-----------------------------------\n";
	std::exit(0);
}

void oneKernel(const string& prog)
{
	std::cerr << "Only one kernel can be specified" << std::endl;
	exitUsage(prog);
}

CommandArgs parse_command_line(int argc, char** argv)
{
	CommandArgs cmd;
	if (argc > 1) cmd.automatic = true;
	for(int a = 1; a < argc; a++)
	{
		auto consume_arg = [&]() -> std::string
		{
			if(a + 1 < argc)
			{
				std::string tmp = argv[a + 1];
				char c = '\0';
				argv[a] = &c; argv[++a] = &c;
				return tmp;
			}
			printf("Warning: missing argument for command '%s'\n", argv[a]);
			return "";
		};
		auto retreive_multiple_gaussian_kernel = [&]() -> void
		{
			if (cmd.customKernel)
				oneKernel(argv[0]);
			if(a + 4 < argc)
			{
				cmd.cutoff = stof(argv[a+1]);
				cmd.kernel.type = Kernel::GAUSSIAN_MULTIPLE;
				cmd.kernel.beginningSigma = stof(argv[a+2]);
				cmd.kernel.factor = stof(argv[a+3]);
				cmd.kernel.nbOfComputedSigmas = unsigned(stoi(argv[a+4]));
				cmd.customKernel = true;
				return;
			}
			exitUsage(argv[0]);
		};
		auto retreive_gaussian_kernel = [&]() -> void
		{
			if (cmd.customKernel)
				oneKernel(argv[0]);
			if(a + 2 < argc)
			{
				cmd.cutoff = stof(argv[a+1]);
				cmd.kernel.type = Kernel::GAUSSIAN;
				cmd.kernel.sigma = stof(argv[a+2]);
				cmd.customKernel = true;
				return;
			}
			exitUsage(argv[0]);
		};
		auto retreive_rational_kernel = [&]() -> void
		{
			if (cmd.customKernel)
				oneKernel(argv[0]);
			if(a + 3 < argc)
			{
				cmd.cutoff = stof(argv[a+1]);
				cmd.kernel.type = Kernel::SINGULAR;
				cmd.kernel.exponent =  stof(argv[a+2]);
				cmd.kernel.constant_shift =  stof(argv[a+3]);
				cmd.customKernel = true;
				return;
			}
			exitUsage(argv[0]);
		};
		auto retreive_unsigned = [&]() -> unsigned
		{
			if (a + 1 < argc)
			{
				return unsigned(stoi(argv[a+1]));
			}
			exitUsage(argv[0]);
			return 0;
		};
		auto retreive_float = [&]() -> float
		{
			if (a + 1 < argc)
			{
				return (stof(argv[a+1]));
			}
			exitUsage(argv[0]);
			return 0;
		};
		if(strcmp(argv[a], "-h") == 0) exitUsage(argv[0]);
		else if(strcmp(argv[a], "-i") == 0) cmd.input_file = consume_arg();
		else if(strcmp(argv[a], "-o") == 0) cmd.output_folder = consume_arg();
		else if(strcmp(argv[a], "-f") == 0) cmd.output_file = consume_arg();
		else if(strcmp(argv[a], "-c") == 0) cmd.cpu = true;
		else if(strcmp(argv[a], "-d") == 0) cmd.dc = true;
		else if(strcmp(argv[a], "-a") == 0) {cmd.knn = true; cmd.nb_knn = retreive_unsigned();}
		else if(strcmp(argv[a], "-e") == 0) {cmd.ball = true; cmd.ball_radius = retreive_float();}
		else if(strcmp(argv[a], "-k") == 0) retreive_multiple_gaussian_kernel();
		else if(strcmp(argv[a], "-r") == 0) retreive_rational_kernel();
		else if(strcmp(argv[a], "-s") == 0) retreive_gaussian_kernel();
		else if(strcmp(argv[a], "-g") == 0) cmd.fast = false;
		else if(strcmp(argv[a], "-m") == 0) cmd.depth = retreive_unsigned();
		else if(strcmp(argv[a], "-b") == 0) cmd.inputOctree = true;
}
return cmd;
}

int main(int argc, char** argv)
{
	srand(unsigned(time(nullptr)));

	const CommandArgs cmd = parse_command_line(argc, argv);
	if(cmd.automatic && cmd.input_file == "")
	{
		std::cerr << "Error, an input file is needed for automatic mode\n";
		exitUsage(argv[0]);
	}

	QApplication a(argc, argv);
	MainWindow w(cmd);
	w.show();

	return a.exec();
}
