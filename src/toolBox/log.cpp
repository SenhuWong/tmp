#include"log.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"

namespace ToolKit
{
	std::shared_ptr<spdlog::logger> Log::s_ConsoleLogger;
	std::shared_ptr<spdlog::logger> Log::s_FileLogger;

	void Log::Init()
	{
		spdlog::set_pattern("%^[%T] %n: %v%$");
		s_ConsoleLogger = spdlog::stdout_color_mt("console");
		s_ConsoleLogger->set_level(spdlog::level::trace);

		s_FileLogger = spdlog::basic_logger_mt("file", "logs/basic_log.txt");
		s_FileLogger->set_level(spdlog::level::trace);
		
	}
	
	

}