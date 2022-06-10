#pragma once
#include <memory>
#include "edge3d_int.h"
#include "spdlog/spdlog.h"
namespace ToolKit
{
	class Log
	{
	private:
		static std::shared_ptr<spdlog::logger> s_ConsoleLogger;
		static std::shared_ptr<spdlog::logger> s_FileLogger;
	public:
		static void Init();
		
		inline static std::shared_ptr<spdlog::logger>& getConsoleLogger() { return s_ConsoleLogger; }
		inline static std::shared_ptr<spdlog::logger>& getFileLogger() { return s_FileLogger; }

	};
}


#ifdef DEBUG
#define TK_CONSOLE_TRACE(...) ::ToolKit::Log::getConsoleLogger()->trace(__VA_ARGS__)
// Core log macros
#define TK_CONSOLE_TRACE(...)    ::ToolKit::Log::getConsoleLogger()->trace(__VA_ARGS__)
#define TK_CONSOLE_INFO(...)     ::ToolKit::Log::getConsoleLogger()->info(__VA_ARGS__)
#define TK_CONSOLE_WARN(...)     ::ToolKit::Log::getConsoleLogger()->warn(__VA_ARGS__)
#define TK_CONSOLE_ERROR(...)    ::ToolKit::Log::getConsoleLogger()->error(__VA_ARGS__)
#define TK_CONSOLE_CRITICAL(...) ::ToolKit::Log::getConsoleLogger()->critical(__VA_ARGS__)

// Client log macros
#define TK_FILE_TRACE(...)         ::ToolKit::Log::getFileLogger()->trace(__VA_ARGS__)
#define TK_FILE_INFO(...)          ::ToolKit::Log::getFileLogger()->info(__VA_ARGS__)
#define TK_FILE_WARN(...)          ::ToolKit::Log::getFileLogger()->warn(__VA_ARGS__)
#define TK_FILE_ERROR(...)         ::ToolKit::Log::getFileLogger()->error(__VA_ARGS__)
#define TK_FILE_CRITICAL(...)      ::ToolKit::Log::getFileLogger()->critical(__VA_ARGS__)
#else
#define TK_CONSOLE_TRACE(...) 
// Core log macros
#define TK_CONSOLE_TRACE(...)    
#define TK_CONSOLE_INFO(...)     
#define TK_CONSOLE_WARN(...)     
#define TK_CONSOLE_ERROR(...)   
#define TK_CONSOLE_CRITICAL(...) 

// Client log macros
#define TK_FILE_TRACE(...)         
#define TK_FILE_INFO(...)          
#define TK_FILE_WARN(...)          
#define TK_FILE_ERROR(...)         
#define TK_FILE_CRITICAL(...)
#endif // DEBUG

template<typename OStream, int ndim, typename T>
inline OStream& operator<<(OStream& os,GeomElements::vector3d<ndim, T>& vec)
{
	return os << vec.to_string();
}

template<typename OStream, int ndim, typename T>
inline OStream& operator<<(OStream& os,const GeomElements::vector3d<ndim, T>& vec)
{
	return os << vec.to_string();
}

template<typename OStream, int ndim>
inline OStream& operator<<(OStream& os,GeomElements::point3d<ndim>& pnt)
{
	return os << pnt.position().to_string();
}

template<typename OStream, int ndim>
inline OStream& operator<<(OStream& os,const GeomElements::point3d<ndim>& pnt)
{
	return os << pnt.position().to_string();
}

template<typename OStream, int ndim>
inline OStream& operator<<(OStream& os, GeomElements::edge3d<ndim>& edg)
{
	os << " Edge size: " << edg.size() << '\n';
	for (int i = 0; i < edg.size(); i++)
	{
		os << i << "th :";
		os << *(edg.point(i));
	}

	return os;
}

template<typename OStream, int ndim>
inline OStream& operator<<(OStream& os,const GeomElements::edge3d<ndim>& edg)
{
	os << " Edge size: " << edg.size() << '\n';
	for (int i = 0; i < edg.size(); i++)
	{
		os << i << "th :";
		os << *(edg.point(i));
	}

	return os;
}