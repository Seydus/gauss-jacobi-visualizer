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
          50: "#EEF6FC",
          100: "#D9E8F6",
          200: "#BCD7EF",
          500: "#427AB5",
          600: "#356291",
          700: "#284F78"
        }
      },
      boxShadow: {
        soft: "0 18px 45px rgba(66, 122, 181, 0.08)"
      }
    }
  }
};
