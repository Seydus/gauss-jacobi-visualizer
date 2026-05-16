/** @type {import('tailwindcss').Config} */
module.exports = {
  content: ["./app.R"],
  prefix: "tw-",
  corePlugins: {
    preflight: false
  },
  theme: {
    extend: {
      colors: {
        ink: "#1f2937",
        muted: "#64748b",
        brand: {
          50: "#eef6ff",
          100: "#d9ebff",
          200: "#b8d9ff",
          500: "#2c5f93",
          600: "#244f7c",
          700: "#1f4268"
        },
        success: {
          50: "#edf7f1",
          500: "#1f6f54",
          600: "#195d46"
        }
      },
      boxShadow: {
        soft: "0 18px 45px rgba(31, 66, 104, 0.08)"
      }
    }
  }
};
